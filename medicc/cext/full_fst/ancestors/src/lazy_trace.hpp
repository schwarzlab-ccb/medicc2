// REVIEW NOTE: Older downward resolver support. The new entry is am_direct_resolve.
// Recover a single conditional optimum directly from lazy subset residuals.
#include "lazy_engine.hpp"
namespace lazy_cascade {
static int scan(Language &q, const string &s, int cap) {
    if (q.start < 0 || int(s.size()) != q.len)
        return INF;
    int state = q.start, cost = q.offset;
    if (cost + q.h(state) > cap)
        return INF;
    for (char ch : s) {
        int c = ch == 'X' ? 9 : ch - '0';
        auto e = q.step(state, c);
        if (e.to < 0)
            return INF;
        cost += e.w;
        state = e.to;
        if (cost + q.h(state) > cap)
            return INF;
    }
    return min(INF, cost + q.final(state));
}
static string recover(Language &q, const string &input) {
    if (dynamic_cast<Fixed *>(&q))
        return input;
    if (auto *h = dynamic_cast<Half *>(&q)) {
        string next = input;
        for (char &c : next)
            if (c != 'X') {
                c = '0' + wgd_model::double_cn(c - '0');
            }
        return recover(*h->inner, next);
    }
    auto *p = dynamic_cast<Pull *>(&q);
    if (!p)
        throw runtime_error("invalid lazy traceback stage");
    vector<int> path = {p->start}, weights;
    int state = p->start;
    for (char ch : input) {
        auto e = p->step(state, ch == 'X' ? 9 : ch - '0');
        if (e.to < 0)
            throw runtime_error("missing lazy traceback arc");
        weights.push_back(e.w);
        path.push_back(e.to);
        state = e.to;
    }
    int tp = -1, tq = -1, residual = 0, best = INF;
    const auto &last = p->states[state].subset;
    for (size_t k = 0; k < last.size(); k += 3) {
        int cost = last[k + 2] + p->p->fin[last[k]] + p->inner->final(last[k + 1]);
        if (cost < best) {
            best = cost;
            tp = last[k];
            tq = last[k + 1];
            residual = last[k + 2];
        }
    }
    if (tp < 0)
        throw runtime_error("no final lazy traceback state");
    string next(input.size(), '0');
    for (int i = int(input.size()) - 1; i >= 0; i--) {
        int c = input[i] == 'X' ? 9 : input[i] - '0';
        bool found = false;
        const auto &prev = p->states[path[i]].subset;
        for (size_t k = 0; k < prev.size() && !found; k += 3)
            for (const auto &r : p->p->arcs[prev[k]][c])
                if (r.to == tp) {
                    auto e = p->inner->step(prev[k + 1], r.out);
                    if (e.to == tq && prev[k + 2] + r.w + e.w == weights[i] + residual) {
                        next[i] = r.out == 9 ? 'X' : char('0' + r.out);
                        tp = prev[k];
                        tq = prev[k + 1];
                        residual = prev[k + 2];
                        found = true;
                        break;
                    }
                }
        if (!found)
            throw runtime_error("no tight predecessor in lazy traceback");
    }
    return recover(*p->inner, next);
}
} // namespace lazy_cascade
extern "C" const char *am_lazy_resolve(void *h, const char *parent, void *c, int upper, int flags,
                                       int limit, int *score) {
    try {
        using namespace lazy_cascade;
        auto *m = (AM *)h;
        auto &child = *(DFA *)c;
        string x = parent;
        if (limit != 0)
            throw runtime_error("production build requires state_limit=0 (unlimited)");
        Ctx ctx{0, 0, 0, limit, !(flags & 2)};
        ctx.owner = m;
        if ((flags & 4) && upper < INF)
            ctx.budget = upper - child.offset;
        Cascade cascade(m, &ctx, child, false, INF);
        int best = upper;
        Language *winner = nullptr;
        for (auto b : cascade.branches) {
            int cap = best >= INF ? INF : best - b.wg;
            int cost = scan(*b.q, x, cap);
            if (cost < INF && cost + b.wg <= best && (!winner || cost + b.wg < best)) {
                winner = b.q;
                best = cost + b.wg;
            }
        }
        if (!winner)
            throw runtime_error("No feasible conditional ancestor");
        resolved = recover(*winner, x);
        *score = best;
        return resolved.c_str();
    } catch (exception &e) {
        ((AM *)h)->core.budget = INF;
        error = e.what();
        return nullptr;
    }
}
