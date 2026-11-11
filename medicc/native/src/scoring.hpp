#pragma once

// Scalar scores. score() splits a pair of CN strings at X, shortens repeated
// columns per chromosome, and for every WGD branch combination sums exact
// per-chromosome scores (fixed(), cached per handle) in bound order. WGDs are
// shared across chromosomes and charged once. directed=true gives the branch
// score parent -> child; otherwise the symmetric pair distance via a common
// ancestor.

// A complete path supplies a feasible upper bound, never an approximate answer.
static double greedy_pair(const PairGraph &graph) {
    // Keep active_g unset: greedy prefixes are not necessarily shortest paths.
    double cost = 0;
    for (int s = graph.exact->Start(); s >= 0;) {
        double final = graph.exact->Final(s).Value();
        if (std::isfinite(final))
            return cost + final <= graph.budget->cap ? cost + final : infinity;
        double best = infinity, weight = 0;
        int next = -1;
        for (fst::ArcIterator<Fst> it(*graph.exact, s); !it.Done(); it.Next()) {
            const auto &a = it.Value();
            double value = a.weight.Value() + (*graph.heuristic)(a.nextstate);
            if (value < best) {
                best = value;
                weight = a.weight.Value();
                next = a.nextstate;
            }
        }
        if (next < 0 || cost + best > graph.budget->cap)
            return infinity;
        cost += weight;
        s = next;
    }
    return infinity;
}
static double fixed(Model &m, const std::string &x, const std::string &y, int k, int j,
                    bool directed, double cap) {
    if (!directed && (x > y || (x == y && k > j)))
        return fixed(m, y, x, j, k, false, cap);
    const std::string key = score_key(x, y, k, j, directed);
    const auto known = m.chromosome_scores.find(key);
    if (known != m.chromosome_scores.end())
        return known->second;
    double value;
    if (!directed) {
        auto graph = make_pair_graph(m, x, y, k, j, cap);
        double feasible = greedy_pair(graph);
        cap = std::min(cap, feasible - 1);
        graph.budget->cap = cap;
        value = std::min(feasible, search(*graph.exact, cap, false, graph.heuristic.get(),
                                          graph.budget.get()).cost);
    } else {
        Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
        auto graph = make_graph(m, a, b, k, j, directed);
        value = search(*graph, cap, false).cost;
    }
    // A bounded failure proves only value > cap, never global infeasibility.
    if (std::isfinite(value))
        cache_score(m, key, value);
    return value;
}
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    // Fixed WGD counts are shared across chromosomes and charged once.
    std::vector<std::pair<std::string, std::string>> segments(1);
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] == 'X') {
            segments.emplace_back();
            continue;
        }
        auto &s = segments.back();
        if (s.first.empty() || s.first.back() != x[i] || s.second.back() != y[i]) {
            s.first += x[i];
            s.second += y[i];
        }
    }
    std::stable_sort(segments.begin(), segments.end(),
                     [](const auto &a, const auto &b) { return a.first.size() < b.first.size(); });
    std::vector<BlockBounds> bounds(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        bounds[i].upper.fill(infinity);
        auto &[a, b] = segments[i];
        for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k)
            for (int j = 0; j <= m.max_wgds; ++j) {
                int branch = 3 * k + j;
                double known = cached_score(m, a, b, k, j, directed);
                if (known >= 0) {
                    bounds[i].lower[branch] = bounds[i].upper[branch] = known;
                    continue;
                }
                double lower = directed ? 0 : string_lower(a, b, branch),
                       upper = witness(a, b, k, j, directed);
                bounds[i].lower[branch] = lower;
                bounds[i].upper[branch] = upper;
                if (lower == upper)
                    cache_score(m, score_key(a, b, k, j, directed), lower);
            }
    }
    double best = infinity;
    auto branches = branch_order(m, bounds, directed, best);
    for (auto [lower, branch] : branches) {
        if (lower >= best)
            continue;
        double bound = lower;
        int k = branch / 3, j = branch % 3;
        for (size_t i = 0; i < segments.size(); ++i) {
            if (bound >= best)
                break;
            double before = bounds[i].lower[branch], other = bound - before;
            // Integral costs: only a strict improvement is needed for a scalar.
            double room = std::min(best - other - 1, bounds[i].upper[branch] - 1);
            double value = fixed(m, segments[i].first, segments[i].second, k, j, directed, room);
            // No cheaper path certifies the feasible history as optimal.
            if (!std::isfinite(value) && room == bounds[i].upper[branch] - 1) {
                value = bounds[i].upper[branch];
                cache_score(m, score_key(segments[i].first, segments[i].second, k, j, directed),
                            value);
            }
            if (!std::isfinite(value)) {
                bound = infinity;
                break;
            }
            bound += value - before;
        }
        best = std::min(best, bound);
    }
    return best;
}
