// Shared AM lifetime, branch union, and C/Python interoperability.
// The main direct engine includes this through older support files; those older
// entry points are retained for validation, not selected by review/model.py.
#ifndef ANCESTRAL_MAX_WGDS
#define ANCESTRAL_MAX_WGDS 3
#endif
#include "budget_2m_core.hpp"
#include <queue>
#include <tuple>
// All APIs operate on whole genomes. X is an explicit tenth symbol.
// AM owns its primitive Parts by value. Side objects point into these members.
// Therefore AM must remain alive until its searches finish. Never delete a Part
// through a borrowed Side pointer.
struct AM {
    Model core;
    Part loh, lg, l, g, w, d, a, ad;
    double visited = 0, peak = 0;
};
static vector<vector<int>> levels(const DFA &a) {
    vector<vector<int>> ls(a.len + 1);
    if (a.start < 0)
        return ls;
    vector<int> depth(a.v.size(), -1);
    depth[a.start] = 0;
    ls[0].push_back(a.start);
    for (int i = 0; i <= a.len; i++)
        for (int s : ls[i])
            for (int c = 0; c < 10; c++)
                if (a.v[s].to[c] >= 0) {
                    int t = a.v[s].to[c];
                    if (i == a.len)
                        throw runtime_error("arc beyond final layer");
                    if (depth[t] < 0) {
                        depth[t] = i + 1;
                        ls[i + 1].push_back(t);
                    } else if (depth[t] != i + 1)
                        throw runtime_error("nonlayered DFA");
                }
    return ls;
}
static vector<int> suffix(const DFA &a) {
    auto ls = levels(a);
    vector<int> h(a.v.size(), INF);
    for (int i = a.len; i >= 0; i--)
        for (int s : ls[i]) {
            h[s] = a.v[s].fin;
            for (int c = 0; c < 10; c++)
                if (a.v[s].to[c] >= 0)
                    h[s] = min(h[s], a.v[s].w[c] + h[a.v[s].to[c]]);
        }
    return h;
}
static Part inverse(Part p) {
    auto q = p;
    for (auto &x : q.arcs)
        for (auto &y : x)
            y.clear();
    for (int s = 0; s < (int)p.arcs.size(); s++)
        for (int c = 0; c < 10; c++)
            for (auto a : p.arcs[s][c])
                q.arcs[s][a.out].push_back({c, a.w, a.to});
    return q;
}
static vector<DFA> branches(AM *m, const DFA &a, bool forward, bool literal, int cap = INF) {
    vector<DFA> out;
    // Reuse backward LG/D suffixes across globally fixed WGD counts.
    // A stage built at the preceding (larger) budget preserves every path
    // required by the next smaller budget; no stricter cutoff is introduced.
    if (!forward && !literal) {
        m->core.budget = cap;
        out.push_back(pull(&m->core, m->l, pull(&m->core, m->g, a)));
        DFA q = pull(&m->core, m->lg, a);
        for (int k = 0; k <= ANCESTRAL_MAX_WGDS; k++) {
            m->core.budget = cap >= INF ? INF : cap - k;
            if (k > 0) {
                // Capped x2 is a deterministic relabel, including identity X.
                DFA half = q;
                for (auto &n : half.v) {
                    Node v;
                    v.fin = n.fin;
                    for (int c = 0; c <= 8; c++) {
                        v.to[c] = n.to[wgd_model::double_cn(c)];
                        v.w[c] = n.w[wgd_model::double_cn(c)];
                    }
                    v.to[9] = n.to[9];
                    v.w[9] = n.w[9];
                    n = v;
                }
                auto ls = levels(half);
                half = minimize(move(half), ls);
                q = pull(&m->core, m->lg, half);
            }
            DFA branch = pull(&m->core, m->loh, q);
            if (branch.start >= 0)
                branch.offset += k;
            out.push_back(move(branch));
        }
        m->core.budget = INF;
        return out;
    }
    auto apply = [&](DFA x, vector<const Part *> chain) {
        if (!forward)
            reverse(chain.begin(), chain.end());
        for (auto p : chain) {
            if (forward) {
                auto q = inverse(*p);
                x = pull(&m->core, q, x);
            } else
                x = pull(&m->core, *p, x);
        }
        return x;
    };
    m->core.budget = cap;
    out.push_back(apply(a, {&m->l, &m->g}));
    if (literal) {
        vector<const Part *> chain = {&m->loh, &m->lg};
        for (int j = 0; j < ANCESTRAL_MAX_WGDS; j++) {
            chain.push_back(&m->w);
            chain.push_back(&m->lg);
        }
        out.push_back(apply(a, chain));
    } else
        for (int k = 0; k <= ANCESTRAL_MAX_WGDS; k++) {
            vector<const Part *> chain = {&m->loh, &m->lg};
            for (int j = 0; j < k; j++) {
                chain.push_back(&m->d);
                chain.push_back(&m->lg);
            }
            m->core.budget = cap >= INF ? INF : cap - k;
            auto b = apply(a, chain);
            if (b.start >= 0)
                b.offset += k;
            out.push_back(move(b));
        }
    m->core.budget = INF;
    return out;
}
// Weighted union, followed by exact acyclic minimization. No string enumeration.
static DFA unite(const vector<DFA> &ds) {
    auto d=library_language::united(ds);
    auto layers=levels(d);
    return minimize(move(d),layers);
}
// A* explores f <= optimum, then keeps every tight path, including all ties.
static DFA optimal_product(AM *m, const DFA &a, const DFA &b, int cap = INF) {
    DFA out;
    out.len = a.len;
    if (a.start < 0 || b.start < 0)
        return out;
    if (a.len != b.len)
        throw runtime_error("length mismatch");
    auto ha = suffix(a), hb = suffix(b);
    int offset = a.offset + b.offset;
    if (offset + ha[a.start] + hb[b.start] > cap)
        return out;
    struct S {
        int a, b, lev, g = INF, tail = INF;
        Node q;
    };
    vector<S> ss;
    unordered_map<uint64_t, int> ids;
    auto insert = [&](int x, int y, int lev) {
        uint64_t key = (uint64_t(uint32_t(x)) << 32) | uint32_t(y);
        auto it = ids.find(key);
        if (it != ids.end())
            return it->second;
        int id = ss.size();
        ids[key] = id;
        ss.push_back({x, y, lev});
        return id;
    };
    using Item = tuple<int, int, int>;
    priority_queue<Item, vector<Item>, greater<Item>> pq;
    int st = insert(a.start, b.start, 0);
    ss[st].g = 0;
    pq.push({ha[a.start] + hb[b.start], 0, st});
    int best = cap >= INF ? INF : cap - offset;
    while (!pq.empty()) {
        auto [f, g, s] = pq.top();
        pq.pop();
        if (f > best)
            break;
        if (g != ss[s].g)
            continue;
        int x = ss[s].a, y = ss[s].b, lev = ss[s].lev;
        m->visited++;
        auto &an = a.v[x];
        auto &bn = b.v[y];
        if (lev == a.len) {
            int fin = an.fin + bn.fin;
            if (fin < INF) {
                ss[s].q.fin = fin;
                best = min(best, g + fin);
            }
            continue;
        }
        for (int c = 0; c < 10; c++)
            if (an.to[c] >= 0 && bn.to[c] >= 0) {
                int x2 = an.to[c], y2 = bn.to[c], w = an.w[c] + bn.w[c];
                int h = ha[x2] + hb[y2];
                if (g + w + h > best)
                    continue;
                int t = insert(x2, y2, lev + 1);
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
    for (int s = 0; s < (int)ss.size(); s++)
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
    for (int s = 0; s < (int)ss.size(); s++) {
        if (ss[s].g + ss[s].tail != best)
            continue;
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
// C ABI ownership contract:
//   am_load -> new AM, released exactly once by am_free.
//   am_string/am_import/am_nfa_join -> new DFA, released once by am_drop.
//   am_arc/am_final/am_size borrow DFA pointers and do not retain them.
//   am_error/traceback strings are borrowed thread-local buffers, not owned results.
// make_unique protects construction against exceptions. release() transfers the
// raw pointer to Python; it deliberately stops unique_ptr from deleting it here.
// static_cast<AM*>(void_handle) recovers the type; it does not copy/allocate.
// Deleting a DFA also frees its vector members via their C++ destructors.
extern "C" {
const char *am_error() {
    return error.c_str();
}
void *am_load(const char *file) {
    try {
        auto m = make_unique<AM>();
        ifstream f(file);
        m->loh = readpart(f, "LOH");
        m->lg = readpart(f, "LG");
        m->l = readpart(f, "L_LOH");
        m->g = readpart(f, "G");
        m->w = readpart(f, "W1");
        m->a = readpart(f, "A");
        m->ad = readpart(f, "AD");
        m->d.start = 0;
        m->d.fin = {0};
        m->d.arcs.resize(1);
        for (int c = 0; c <= 8; c++)
            m->d.arcs[0][c].push_back({wgd_model::double_cn(c), 0, 0});
        m->d.arcs[0][9].push_back({9, 0, 0});
        return m.release();
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}
void am_free(void *m) {
    delete static_cast<AM *>(m);
}
void am_drop(void *d) {
    delete static_cast<DFA *>(d);
}
void *am_string(const char *s) {
    try {
        return new DFA(linear(s));
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}
void *am_join(void *h, void *l, void *r, int literal, int upper) {
    try {
        auto *m = (AM *)h;
        int ca = upper >= INF ? INF : upper - ((DFA *)r)->offset,
            cb = upper >= INF ? INF : upper - ((DFA *)l)->offset;
        auto as = branches(m, *(DFA *)l, false, literal, ca),
             bs = branches(m, *(DFA *)r, false, literal, cb);
        vector<DFA> kept;
        int best = upper;
        for (auto &a : as)
            for (auto &b : bs) {
                auto d = optimal_product(m, a, b, best);
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
        error = e.what();
        return nullptr;
    }
}
void *am_align(void *h, void *p, void *c, int literal, int upper) {
    try {
        auto *m = (AM *)h;
        int cap = upper >= INF ? INF : upper - ((DFA *)c)->offset;
        auto as = branches(m, *(DFA *)p, true, literal, cap);
        vector<DFA> kept;
        int best = upper;
        for (auto &a : as) {
            auto d = optimal_product(m, a, *(DFA *)c, best);
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
        error = e.what();
        return nullptr;
    }
}
// rows is a borrowed contiguous array of 21 integers per state:
// [final_weight, next_for_0, weight_for_0, ..., next_for_X, weight_for_X].
// am_import copies these numbers into owned vectors before returning, so Python
// may release the input array after the call. It does not save rows itself.
void *am_import(int len, int n, int start, int offset, const int *rows) {
    try {
        if (len < 0 || n < 0 || start < -1 || start >= n || offset < 0 || offset >= INF)
            throw runtime_error("invalid imported DFA header");
        DFA d;
        d.len = len;
        d.start = start;
        d.offset = offset;
        d.v.resize(n);
        for (int s = 0; s < n; s++) {
            auto &q = d.v[s];
            q.fin = rows[21 * s];
            if (q.fin < 0 || q.fin > INF)
                throw runtime_error("invalid final weight");
            for (int c = 0; c < 10; c++) {
                q.to[c] = rows[21 * s + 1 + 2 * c];
                q.w[c] = rows[21 * s + 2 + 2 * c];
                if (q.to[c] < -1 || q.to[c] >= n || q.w[c] < 0 || q.w[c] >= INF)
                    throw runtime_error("invalid imported arc");
            }
        }
        auto ls = levels(d);
        for (int i = 0; i < len; i++)
            for (int s : ls[i])
                if (d.v[s].fin < INF)
                    throw runtime_error("variable-length acceptor");
        return new DFA(minimize(move(d), ls));
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}
int am_size(void *d) {
    return ((DFA *)d)->v.size();
}
int am_start(void *d) {
    return ((DFA *)d)->start;
}
int am_offset(void *d) {
    return ((DFA *)d)->offset;
}
int am_arc(void *d, int s, int c, int *w) {
    auto &a = ((DFA *)d)->v.at(s);
    *w = a.w.at(c);
    return a.to.at(c);
}
int am_final(void *d, int s) {
    return ((DFA *)d)->v.at(s).fin;
}
void am_stats(void *h, double *out) {
    auto *m = (AM *)h;
    out[0] = m->core.stats[1];
    out[1] = m->core.stats[2];
    out[2] = m->visited;
    out[3] = m->peak;
}
}
