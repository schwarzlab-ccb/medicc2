#pragma once

// Minimum-cost event histories for one parent/child branch (used by --events).
// history() recovers the optimal WGD count, then traces one optimal path through
// the composed stages A (D A)^k and splits each A into its loss and gain stage.
// The recorded composition tuples recover intermediate CNs on that path; OpenFST
// owns every product state and no unrestricted P or history transducer is built.
// Output: the score, then one "type<TAB>before<TAB>after" line per stage.
struct StageTrace {
    double cost;
    std::vector<std::string> words;
};
static StageTrace trace_stages(const std::string &x, const std::string &y,
                               const std::vector<const Vec *> &stages, double cap) {
    Ptr graph(new Vec(linear(y)));
    std::vector<StateTable *> tables(stages.size());
    for (int i = int(stages.size()) - 1; i >= 0; --i)
        graph = tracked_compose(*stages[i], *graph, &tables[i]);
    auto parent = linear(x);
    StateTable *outer;
    graph = tracked_compose(parent, *graph, &outer);
    std::vector<int> ids;
    auto result = search(*graph, cap, true, nullptr, nullptr, &ids);
    if (!std::isfinite(result.cost))
        throw std::runtime_error("No optimal history within score bound");
    std::vector<W> backward;
    fst::ShortestDistance(result.tight, &backward, true);
    StageTrace trace{result.cost, std::vector<std::string>(stages.size() + 1)};
    trace.words[0] = x;
    int state = result.tight.Start();
    size_t column = 0;
    while (result.tight.Final(state) == W::Zero()) {
        bool found = false;
        for (fst::ArcIterator<Vec> it(result.tight, state); !it.Done(); it.Next()) {
            auto arc = it.Value();
            if (!std::isfinite(distance_at(backward, arc.nextstate)))
                continue;
            int source = outer->Tuple(ids[state]).StateId2();
            int target = outer->Tuple(ids[arc.nextstate]).StateId2();
            // Decode a synchronous composed arc by a tiny label DP. These are
            // labels on already selected native states, not new product states.
            std::vector<std::array<double, 11>> costs(stages.size() + 1);
            std::vector<std::array<int, 11>> previous(stages.size() + 1);
            for (auto &row : costs)
                row.fill(infinity);
            costs[0][label(x.at(column))] = 0;
            for (size_t i = 0; i < stages.size(); ++i) {
                auto a = tables[i]->Tuple(source), b = tables[i]->Tuple(target);
                for (fst::ArcIterator<Vec> ai(*stages[i], a.StateId1()); !ai.Done(); ai.Next()) {
                    const auto &v = ai.Value();
                    if (v.nextstate != b.StateId1())
                        continue;
                    double cost = costs[i][v.ilabel] + v.weight.Value();
                    if (cost < costs[i + 1][v.olabel]) {
                        costs[i + 1][v.olabel] = cost;
                        previous[i + 1][v.olabel] = v.ilabel;
                    }
                }
                source = a.StateId2();
                target = b.StateId2();
            }
            int value = label(y.at(column));
            if (!std::isfinite(costs.back()[value]))
                throw std::runtime_error("Cannot decode native history arc");
            for (int i = int(stages.size()); i > 0; --i) {
                trace.words[i] += symbol(value);
                value = previous[i][value];
            }
            state = arc.nextstate;
            ++column;
            found = true;
            break;
        }
        if (!found)
            throw std::runtime_error("Broken optimal history path");
    }
    if (column != x.size())
        throw std::runtime_error("History length mismatch");
    return trace;
}
static std::vector<std::string> chromosome_words(const std::string &s) {
    std::vector<std::string> out(1);
    for (char c : s) {
        if (c == 'X')
            out.emplace_back();
        else
            out.back() += c;
    }
    return out;
}
static std::string history(Model &m, const std::string &x, const std::string &y) {
    compatible(x, y);
    double optimum = score(m, x, y, true);
    if (!std::isfinite(optimum))
        throw std::runtime_error("No finite directed history");
    auto xs = chromosome_words(x), ys = chromosome_words(y);
    int selected = -1;
    std::vector<double> costs(xs.size());
    for (int k = 0; k <= m.max_wgds && selected < 0; ++k) {
        double total = k;
        for (size_t c = 0; c < xs.size() && total <= optimum; ++c) {
            std::string a, b;
            for (size_t i = 0; i < xs[c].size(); ++i)
                if (i == 0 || xs[c][i] != xs[c][i - 1] || ys[c][i] != ys[c][i - 1]) {
                    a += xs[c][i];
                    b += ys[c][i];
                }
            costs[c] = fixed(m, a, b, 0, k, true, optimum - total);
            total += costs[c];
        }
        if (total == optimum)
            selected = k;
    }
    if (selected < 0)
        throw std::runtime_error("Cannot recover optimal WGD count");
    std::vector<const Vec *> stages{&m.a};
    for (int k = 0; k < selected; ++k) {
        stages.push_back(&m.doubling);
        stages.push_back(&m.a);
    }
    std::vector<std::string> full(3 * selected + 3);
    for (size_t c = 0; c < xs.size(); ++c) {
        auto trace = trace_stages(xs[c], ys[c], stages, costs[c]);
        if (trace.cost != costs[c])
            throw std::runtime_error("History and branch score differ");
        std::vector<std::string> words{xs[c]};
        for (size_t i = 0; i < stages.size(); ++i) {
            if (stages[i] == &m.a) {
                auto local = trace_stages(trace.words[i], trace.words[i + 1], {&m.loss, &m.gain},
                                          changes(trace.words[i], trace.words[i + 1]));
                words.push_back(local.words[1]);
                words.push_back(local.words[2]);
            } else
                words.push_back(trace.words[i + 1]);
        }
        for (size_t i = 0; i < full.size(); ++i) {
            if (c)
                full[i] += 'X';
            full[i] += words[i];
        }
    }
    std::string text = std::to_string(int(optimum)) + "\n";
    for (size_t i = 1; i < full.size(); ++i) {
        const char *kind = i % 3 == 1 ? "loss" : i % 3 == 2 ? "gain" : "wgd";
        text += std::string(kind) + "\t" + full[i - 1] + "\t" + full[i] + "\n";
    }
    return text;
}
