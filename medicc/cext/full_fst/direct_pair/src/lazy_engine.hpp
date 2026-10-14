// REVIEW NOTE: Older lazy weighted-subset composition support; not the primary review path.
// Lazy weighted-subset cascade: no complete intermediate ancestor DFA is built.
// Reuse the established adapters and conditional down-pass; replace sibling join.
#include "shared_trace.hpp"
namespace lazy_cascade {
struct Ctx {
    size_t nodes = 0, arcs = 0, entries = 0;
    int limit;
    bool heuristic;
    int budget = INF;
    AM *owner = nullptr;
    ~Ctx() {
        if (owner) {
            owner->core.stats[5] = max(owner->core.stats[5], double(nodes));
            owner->core.stats[6] = max(owner->core.stats[6], double(entries));
            owner->core.stats[7] += arcs;
        }
    }
};
struct Edge {
    int to = -1, w = 0;
};
struct Language {
    int start = -1, len = 0, offset = 0;
    virtual ~Language() = default;
    virtual Edge step(int, int) = 0;
    virtual int final(int) = 0;
    virtual int h(int) = 0;
    virtual int even_h(int s) {
        return h(s);
    }
};
static vector<int> restricted_potential(const Part &p, const DFA &a, bool even) {
    auto ls = levels(a);
    int np = p.fin.size();
    vector<int> h(a.v.size() * np, INF);
    for (int lev = a.len; lev >= 0; lev--)
        for (int s : ls[lev])
            for (int ps = 0; ps < np; ps++) {
                int best = min(INF, p.fin[ps] + a.v[s].fin);
                for (int c = 0; c < 10; c++)
                    if (!even || c == 9 || c % 2 == 0)
                        for (auto &r : p.arcs[ps][c]) {
                            int t = a.v[s].to[r.out];
                            if (t >= 0)
                                best = min(best, r.w + a.v[s].w[r.out] + h[size_t(t) * np + r.to]);
                        }
                h[size_t(s) * np + ps] = min(best, INF);
            }
    return h;
}
struct Fixed : Language {
    DFA owned;
    const DFA *a;
    vector<int> hs, es;
    explicit Fixed(const DFA &d) : a(&d) {
        start = d.start;
        len = d.len;
        offset = d.offset;
        hs = suffix(d);
    }
    explicit Fixed(DFA &&d) : owned(move(d)), a(&owned) {
        start = a->start;
        len = a->len;
        offset = a->offset;
        hs = suffix(*a);
    }
    Edge step(int s, int c) override {
        return {a->v[s].to[c], a->v[s].w[c]};
    }
    int final(int s) override {
        return a->v[s].fin;
    }
    int h(int s) override {
        return hs[s];
    }
    int even_h(int s) override {
        if (es.empty()) {
            auto ls = levels(*a);
            es.assign(a->v.size(), INF);
            for (int lev = len; lev >= 0; lev--)
                for (int x : ls[lev]) {
                    es[x] = a->v[x].fin;
                    for (int c = 0; c < 10; c++)
                        if (c == 9 || c % 2 == 0) {
                            auto e = step(x, c);
                            if (e.to >= 0)
                                es[x] = min(es[x], e.w + es[e.to]);
                        }
                }
        }
        return es[s];
    }
};
struct Pull : Language {
    struct State {
        vector<int> subset;
        array<int, 10> to, w;
        int lev, fin = -2, bound = -2, ebound = -2;
        State(vector<int> v, int l) : subset(move(v)), lev(l) {
            to.fill(-2);
            w.fill(0);
        }
    };
    Ctx *ctx;
    const Part *p;
    Language *inner;
    vector<State> states;
    unordered_map<vector<int>, int, VH> ids;
    vector<int> potential, even_potential;
    int np;
    Pull(Ctx *c, const Part &part, Language *q, bool bottom = false)
        : ctx(c), p(&part), inner(q), np(part.fin.size()) {
        len = q->len;
        offset = q->offset;
        if (q->start < 0)
            return;
        if (bottom && ctx->heuristic) {
            auto *f = dynamic_cast<Fixed *>(q);
            if (!f)
                throw runtime_error("invalid bottom stage");
            potential = restricted_potential(part, *f->a, false);
            even_potential = restricted_potential(part, *f->a, true);
        }
        start = intern({p->start, q->start, 0}, 0);
    }
    int intern(vector<int> sub, int lev) {
        auto key = sub;
        key.push_back(-lev - 1);
        auto it = ids.find(key);
        if (it != ids.end())
            return it->second;
        ++ctx->nodes;
        ctx->entries += sub.size() / 3;
        int id = states.size();
        states.emplace_back(move(sub), lev);
        ids.emplace(move(key), id);
        return id;
    }
    Edge step(int s, int c) override {
        int t = states[s].to[c];
        if (t != -2)
            return {t, states[s].w[c]};
        if (states[s].lev == len) {
            states[s].to[c] = -1;
            return {};
        }
        vector<pair<uint64_t, int>> cand;
        // Inner calls cannot mutate this outer stage. Do not retain refs across intern().
        const auto &sub = states[s].subset;
        for (size_t k = 0; k < sub.size(); k += 3)
            for (auto &r : p->arcs[sub[k]][c]) {
                auto e = inner->step(sub[k + 1], r.out);
                if (e.to < 0)
                    continue;
                int value = sub[k + 2] + r.w + e.w;
                // Fixed call-wide budget only. Never use the prefix that first visits a state.
                // Preserve every stage path whose total residual cost is at most budget.
                if (ctx->budget < INF) {
                    int tail =
                        potential.empty() ? inner->h(e.to) : potential[size_t(e.to) * np + r.to];
                    if (value + tail > ctx->budget)
                        continue;
                }
                uint64_t key = (uint64_t(uint32_t(r.to)) << 32) | uint32_t(e.to);
                cand.push_back({key, value});
            }
        if (cand.empty()) {
            states[s].to[c] = -1;
            return {};
        }
        sort(cand.begin(), cand.end());
        vector<int> out;
        int low = INF;
        for (size_t k = 0; k < cand.size();) {
            size_t end = k + 1;
            int val = cand[k].second;
            while (end < cand.size() && cand[end].first == cand[k].first) {
                val = min(val, cand[end].second);
                end++;
            }
            out.push_back(int(cand[k].first >> 32));
            out.push_back(int(uint32_t(cand[k].first)));
            out.push_back(val);
            low = min(low, val);
            k = end;
        }
        for (size_t k = 2; k < out.size(); k += 3)
            out[k] -= low;
        int lev = states[s].lev;
        t = intern(move(out), lev + 1);
        states[s].to[c] = t;
        states[s].w[c] = low;
        ctx->arcs++;
        return {t, low};
    }
    int final(int s) override {
        if (states[s].lev != len)
            return INF;
        if (states[s].fin != -2)
            return states[s].fin;
        int best = INF;
        auto &sub = states[s].subset;
        for (size_t k = 0; k < sub.size(); k += 3)
            best = min(best, sub[k + 2] + p->fin[sub[k]] + inner->final(sub[k + 1]));
        return states[s].fin = min(best, INF);
    }
    int h(int s) override {
        if (!ctx->heuristic)
            return 0;
        if (states[s].bound != -2)
            return states[s].bound;
        int best = INF;
        auto &sub = states[s].subset;
        for (size_t k = 0; k < sub.size(); k += 3) {
            int v = potential.empty() ? inner->h(sub[k + 1])
                                      : potential[size_t(sub[k + 1]) * np + sub[k]];
            best = min(best, sub[k + 2] + v);
        }
        return states[s].bound = min(best, INF);
    }
    int even_h(int s) override {
        if (!ctx->heuristic || even_potential.empty())
            return h(s);
        if (states[s].ebound != -2)
            return states[s].ebound;
        int best = INF;
        auto &sub = states[s].subset;
        for (size_t k = 0; k < sub.size(); k += 3)
            best = min(best, sub[k + 2] + even_potential[size_t(sub[k + 1]) * np + sub[k]]);
        return states[s].ebound = min(best, INF);
    }
};
struct Half : Language {
    Language *inner;
    explicit Half(Language *q) : inner(q) {
        start = q->start;
        len = q->len;
        offset = q->offset;
    }
    Edge step(int s, int c) override {
        return c == 9 ? inner->step(s, 9) : (c <= 4 ? inner->step(s, 2 * c) : Edge{});
    }
    int final(int s) override {
        return inner->final(s);
    }
    int h(int s) override {
        return inner->even_h(s);
    }
};
struct Branch {
    Language *q;
    int wg;
};
struct Cascade {
    vector<unique_ptr<Language>> owned;
    vector<Branch> branches;
    template <class T, class... Args> T *add(Args &&...args) {
        auto v = make_unique<T>(std::forward<Args>(args)...);
        auto *p = v.get();
        owned.push_back(move(v));
        return p;
    }
    Cascade(AM *m, Ctx *ctx, const DFA &d, bool seeded, int cap) {
        auto *base = add<Fixed>(d);
        Language *g;
        Language *q;
        if (seeded) {
            m->core.budget = cap;
            g = add<Fixed>(pull(&m->core, m->g, d));
            q = add<Fixed>(pull(&m->core, m->lg, d));
            m->core.budget = INF;
        } else {
            g = add<Pull>(ctx, m->g, base, true);
            q = add<Pull>(ctx, m->lg, base, true);
        }
        branches.push_back({add<Pull>(ctx, m->l, g), 0});
        for (int k = 0; k <= ANCESTRAL_MAX_WGDS; k++) {
            if (k > 0)
                q = add<Pull>(ctx, m->lg, add<Half>(q));
            branches.push_back({add<Pull>(ctx, m->loh, q), k});
        }
    }
};
static vector<int> layout(const DFA &d) {
    auto ls = levels(d);
    vector<int> mask(d.len, 0);
    for (int i = 0; i < d.len; i++)
        for (int s : ls[i])
            for (int c = 0; c < 10; c++)
                if (d.v[s].to[c] >= 0)
                    mask[i] |= 1 << c;
    return mask;
}
static DFA product(AM *m, Ctx *ctx, Branch aa, Branch bb, const vector<int> &mask, int cap) {
    auto &a = *aa.q;
    auto &b = *bb.q;
    DFA out;
    out.len = a.len;
    if (a.start < 0 || b.start < 0)
        return out;
    int offset = a.offset + b.offset + aa.wg + bb.wg;
    int h0 = a.h(a.start) + b.h(b.start);
    if (h0 >= INF || offset + h0 > cap)
        return out;
    struct S {
        int a, b, lev, g = INF, tail = INF;
        Node q;
    };
    vector<S> ss;
    unordered_map<uint64_t, int> ids;
    auto insert = [&](int x, int y, int lev) {
        uint64_t k = (uint64_t(uint32_t(x)) << 32) | uint32_t(y);
        auto it = ids.find(k);
        if (it != ids.end())
            return it->second;
        int id = ss.size();
        ss.push_back({x, y, lev});
        ids.emplace(k, id);
        return id;
    };
    using Item = tuple<int, int, int>;
    priority_queue<Item, vector<Item>, greater<Item>> pq;
    int st = insert(a.start, b.start, 0);
    ss[st].g = 0;
    pq.push({h0, 0, st});
    int best = cap >= INF ? INF : cap - offset;
    while (!pq.empty()) {
        auto [f, g, s] = pq.top();
        pq.pop();
        if (f > best)
            break;
        if (g != ss[s].g)
            continue;
        m->visited++;
        int x = ss[s].a, y = ss[s].b, lev = ss[s].lev;
        if (lev == a.len) {
            int fin = a.final(x) + b.final(y);
            if (fin < INF) {
                ss[s].q.fin = fin;
                best = min(best, g + fin);
            }
            continue;
        }
        for (int c = 0; c < 10; c++)
            if (mask[lev] & (1 << c)) {
                auto e = a.step(x, c);
                if (e.to < 0)
                    continue;
                int ha = a.h(e.to);
                if (g + e.w + ha > best)
                    continue;
                auto v = b.step(y, c);
                if (v.to < 0)
                    continue;
                int w = e.w + v.w, h = ha + b.h(v.to);
                if (h >= INF || g + w + h > best)
                    continue;
                int t = insert(e.to, v.to, lev + 1);
                ss[s].q.to[c] = t;
                ss[s].q.w[c] = w;
                if (g + w < ss[t].g) {
                    ss[t].g = g + w;
                    pq.push({g + w + h, g + w, t});
                }
            }
    }
    m->peak = max(m->peak, double(ss.size()));
    if (best >= INF)
        return out;
    vector<vector<int>> ls(a.len + 1);
    for (int s = 0; s < int(ss.size()); s++)
        ls[ss[s].lev].push_back(s);
    for (int i = a.len; i >= 0; i--)
        for (int s : ls[i]) {
            ss[s].tail = ss[s].q.fin;
            for (int c = 0; c < 10; c++)
                if (ss[s].q.to[c] >= 0)
                    ss[s].tail = min(ss[s].tail, ss[s].q.w[c] + ss[ss[s].q.to[c]].tail);
        }
    if (ss[st].tail > best)
        return out;
    out.start = st;
    out.offset = offset;
    out.v.resize(ss.size());
    for (int s = 0; s < int(ss.size()); s++)
        if (ss[s].g + ss[s].tail == best) {
            if (ss[s].g + ss[s].q.fin == best)
                out.v[s].fin = ss[s].q.fin;
            for (int c = 0; c < 10; c++) {
                int t = ss[s].q.to[c];
                if (t >= 0 && ss[s].g + ss[s].q.w[c] + ss[t].tail == best) {
                    out.v[s].to[c] = t;
                    out.v[s].w[c] = ss[s].q.w[c];
                }
            }
        }
    return minimize(move(out), ls);
}
} // namespace lazy_cascade
extern "C" void *am_lazy_join(void *h, void *l, void *r, int upper, int flags, int limit) {
    try {
        using namespace lazy_cascade;
        auto *m = (AM *)h;
        auto &a = *(DFA *)l;
        auto &b = *(DFA *)r;
        if (a.len != b.len)
            throw runtime_error("length mismatch");
        if (limit != 0)
            throw runtime_error("production build requires state_limit=0 (unlimited)");
        if (a.start < 0 || b.start < 0) {
            auto *d = new DFA();
            d->len = a.len;
            return d;
        }
        Ctx ctx{0, 0, 0, limit, !(flags & 2)};
        ctx.owner = m;
        if ((flags & 4) && upper < INF)
            ctx.budget = upper - a.offset - b.offset;
        int ca = upper >= INF ? INF : upper - b.offset, cb = upper >= INF ? INF : upper - a.offset;
        Cascade left(m, &ctx, a, flags & 1, ca), right(m, &ctx, b, flags & 1, cb);
        auto lm = layout(a), rm = layout(b);
        vector<int> mask(a.len);
        for (int i = 0; i < a.len; i++) {
            if ((lm[i] & 511) && (rm[i] & 511))
                mask[i] = 510 | ((lm[i] & 1) && (rm[i] & 1) ? 1 : 0);
            if ((lm[i] & 512) && (rm[i] & 512))
                mask[i] |= 512;
        }
        vector<array<int, 3>> order;
        for (int i = 0; i < int(left.branches.size()); i++)
            for (int j = 0; j < int(right.branches.size()); j++) {
                auto x = left.branches[i], y = right.branches[j];
                if (x.q->start < 0 || y.q->start < 0)
                    continue;
                int lb = x.q->offset + y.q->offset + x.wg + y.wg + x.q->h(x.q->start) +
                         y.q->h(y.q->start);
                order.push_back({lb, i, j});
            }
        sort(order.begin(), order.end());
        int best = upper;
        vector<DFA> kept;
        for (auto item : order) {
            if (item[0] > best)
                continue;
            auto d = product(m, &ctx, left.branches[item[1]], right.branches[item[2]], mask, best);
            if (d.start < 0)
                continue;
            if (d.offset < best) {
                best = d.offset;
                kept.clear();
            }
            if (d.offset == best)
                kept.push_back(move(d));
        }
        return new DFA(unite(kept));
    } catch (exception &e) {
        ((AM *)h)->core.budget = INF;
        error = e.what();
        return nullptr;
    }
}
extern "C" void am_lazy_stats(void *h, double *out) {
    auto *m = (AM *)h;
    out[0] = m->core.stats[5];
    out[1] = m->core.stats[6];
    out[2] = m->core.stats[7];
}
