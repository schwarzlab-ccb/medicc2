// REVIEW NOTE: Shared traceback utilities included by the native compilation unit.
// Keep the retained child language as the terminal constraint during down-pass.
#include "shared_engine.hpp"
static thread_local string resolved;
static int string_cost(const DFA &d, const string &s) {
    if (d.start < 0 || int(s.size()) != d.len)
        return INF;
    int state = d.start, cost = d.offset;
    for (char ch : s) {
        int c = ch == 'X' ? 9 : ch - '0';
        if (c < 0 || c > 9 || d.v[state].to[c] < 0)
            return INF;
        cost += d.v[state].w[c];
        state = d.v[state].to[c];
    }
    return min(INF, cost + d.v[state].fin);
}
static string trace_primitive(const Part &p, const string &input, const DFA &next, int cap) {
    struct R {
        int ps, ds, g, prev, label;
    };
    vector<R> records = {{p.start, next.start, next.offset, -1, -1}};
    vector<int> current = {0};
    auto pot = pull_potential(p, next);
    int np = p.fin.size();
    for (char ch : input) {
        int in = ch == 'X' ? 9 : ch - '0';
        unordered_map<uint64_t, int> memo;
        vector<int> layer;
        for (int id : current) {
            auto r = records[id];
            for (auto &a : p.arcs[r.ps][in]) {
                int t = next.v[r.ds].to[a.out];
                if (t < 0)
                    continue;
                int cost = r.g + a.w + next.v[r.ds].w[a.out];
                if (cost + pot[size_t(t) * np + a.to] > cap)
                    continue;
                uint64_t key = (uint64_t(uint32_t(a.to)) << 32) | uint32_t(t);
                auto it = memo.find(key);
                if (it == memo.end()) {
                    int i = records.size();
                    records.push_back({a.to, t, cost, id, a.out});
                    memo[key] = i;
                    layer.push_back(i);
                } else if (cost < records[it->second].g)
                    records[it->second] = {a.to, t, cost, id, a.out};
            }
        }
        current = move(layer);
    }
    int best = INF, id = -1;
    for (int i : current) {
        auto &r = records[i];
        int cost = r.g + p.fin[r.ps] + next.v[r.ds].fin;
        if (cost < best) {
            best = cost;
            id = i;
        }
    }
    if (id < 0 || best > cap)
        throw runtime_error("primitive traceback has no feasible path");
    string answer;
    for (; records[id].prev >= 0; id = records[id].prev)
        answer.push_back(records[id].label == 9 ? 'X' : char('0' + records[id].label));
    reverse(answer.begin(), answer.end());
    return answer;
}
extern "C" const char *am_resolve(void *h, const char *parent, void *c, int upper, int *score) {
    try {
        auto *m = (AM *)h;
        auto &child = *(DFA *)c;
        string x = parent;
        int best = upper;
        bool found = false;
        for (int k = -1; k <= 3; k++) {
            vector<const Part *> chain;
            if (k < 0)
                chain = {&m->l, &m->g};
            else {
                chain = {&m->loh, &m->lg};
                for (int i = 0; i < k; i++) {
                    chain.push_back(&m->d);
                    chain.push_back(&m->lg);
                }
            }
            int wg = max(0, k);
            m->core.budget = upper >= INF ? INF : best - wg;
            vector<DFA> stages(chain.size() + 1);
            stages.back() = child;
            for (int i = int(chain.size()) - 1; i >= 0; i--)
                stages[i] = pull(&m->core, *chain[i], stages[i + 1]);
            int value = string_cost(stages[0], x);
            if (value >= INF || value + wg > best || (found && value + wg == best))
                continue;
            string s = x;
            for (int i = 0; i < int(chain.size()); i++) {
                int cap = string_cost(stages[i], s);
                s = trace_primitive(*chain[i], s, stages[i + 1], cap);
            }
            best = value + wg;
            resolved = s;
            found = true;
        }
        m->core.budget = INF;
        if (!found)
            throw runtime_error("No feasible conditional ancestor");
        *score = best;
        return resolved.c_str();
    } catch (exception &e) {
        ((AM *)h)->core.budget = INF;
        error = e.what();
        return nullptr;
    }
}
