// Optional event traceback on the SAME (A D)^k A cascade used for scoring.
// The small primitive states along a winning path are enough to recover one
// sequence of intermediate profiles. No eagerly composed genome FST is built.
namespace event_trace {
static vector<int> column(const direct_nfa::Side &side, uint64_t from,
                          uint64_t to, int input, int output) {
    const auto state = side.prepare(from);
    const uint32_t next = to >> 32;
    vector<int> values(side.parts.size()), best;
    int cost = INF;
    auto visit = [&](auto &&visit, int stage, int symbol, int weight) -> void {
        if (weight >= cost) return;
        const auto &part = *side.parts[stage];
        int target = (next / side.radix[stage]) % part.fin.size();
        for (const auto &arc : part.arcs[state.ps[stage]][symbol]) {
            if (arc.to != target) continue;
            values[stage] = arc.out;
            if (stage + 1 == int(side.parts.size())) {
                if (arc.out == output && weight + arc.w < cost) {
                    cost = weight + arc.w;
                    best = values;
                }
            } else visit(visit, stage + 1, arc.out, weight + arc.w);
        }
    };
    visit(visit, 0, input, 0);
    if (best.empty()) throw runtime_error("Cannot recover primitive event transition");
    return best;
}
}

// Returns JSON in a borrowed thread-local buffer (ctypes copies it immediately).
// upper is an exact directed score supplied by the caller; recover an actual
// accepting path within that budget. Prefer fewer WGDs among equal-cost paths.
extern "C" const char *am_pair_history(void *handle, const char *parent,
                                        const char *child, int upper) {
    static thread_local string response;
    try {
        string x = parent, y = child;
        if (x.size() != y.size() || upper < 0 || upper >= INF)
            throw runtime_error("Invalid event traceback bound or length");
        for (size_t i = 0; i < x.size(); ++i) {
            for (char c : {x[i], y[i]})
                if (c != 'X' && (c < '0' || c > '8')) throw runtime_error("Invalid CN symbol");
            if ((x[i] == 'X') != (y[i] == 'X')) throw runtime_error("Separator mismatch");
        }
        DFA target = linear(y);
        auto choices = direct_nfa::branches(static_cast<AM *>(handle), target, 0, false);
        for (auto &choice : choices) {
            auto &side = *choice;
            if (side.wg + side.h(side.start()) > upper) continue;
            vector<uint64_t> path;
            double stats[3] = {};
            auto found = direct_trace::search(side, x, upper, stats, &path);
            if (found.cost >= INF) continue;
            vector<string> stages(side.parts.size(), string(x.size(), 'X'));
            for (size_t i = 0; i < x.size(); ++i) {
                auto labels = event_trace::column(side, path[i], path[i+1],
                    x[i] == 'X' ? 9 : x[i]-'0', y[i] == 'X' ? 9 : y[i]-'0');
                for (size_t j = 0; j < labels.size(); ++j)
                    stages[j][i] = labels[j] == 9 ? 'X' : char('0' + labels[j]);
            }
            ostringstream out;
            out << "{\"cost\":" << found.cost << ",\"wgds\":" << side.wg << ",\"stages\":[";
            for (size_t i = 0; i < stages.size(); ++i) {
                if (i) out << ',';
                out << '"' << stages[i] << '"';
            }
            out << "]}";
            response = out.str();
            return response.c_str();
        }
        throw runtime_error("No capped event history within the directed score");
    } catch (exception &e) { error = e.what(); return nullptr; }
}
