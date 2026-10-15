#if defined(IE_VARIATION_BOUND) || defined(IE_LOCAL_BOUND)
#error "Legacy optional bound tables are for strict doubling; use the capped chromosome bound"
#endif
#include "../../common/src/wgd_model.hpp"
// REVIEW NOTE: Independent exact pair scorer: feasible witnesses and final edge rescoring.
// C ABI ie_load/ie_free own a Model pointer; ie_score borrows it. Its cache
// uses shared_ptr<DFA> (reference-counted ownership). get() on a shared_ptr
// borrows its raw pointer; moving the shared_ptr transfers an owner, not the DFA.
// C strings x/y are borrowed only for the call. Exceptions become a negative
// error code; infinity is a legitimate infeasible distance, not an exception.
// The model holds its WGD-count choice across all chromosomes. Summing independently
// optimized per-chromosome WGD counts would describe a different biological model.
// Integer weighted subset construction and acyclic minimization.
// All new code. Components are exported from the specified MEDICC factory.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef IE_VARIATION_BOUND
#include "generated/variation_bound.hpp"
#endif
#ifdef IE_LOCAL_BOUND
#include "generated/local_bound.hpp"
#endif
#ifndef IE_MAX_WGDS
#define IE_MAX_WGDS 3
#endif
using namespace std;
constexpr int INF = 100000000;
struct Arc {
    int out, w, to;
};
struct Part {
    int start;
    vector<int> fin;
    vector<array<vector<Arc>, 9>> arcs;
};
struct Node {
    array<int, 9> to, w;
    int fin = INF;
    Node() {
        to.fill(-1);
        w.fill(0);
    }
};
struct DFA {
    vector<Node> v;
    int start = -1, offset = 0, len = 0;
};
struct VH {
    size_t operator()(const vector<int> &v) const {
        size_t h = 0x9e3779b97f4a7c15ULL;
        for (int x : v)
            h ^= size_t(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};
struct Model {
    Part loh, lg, now;
    int kind = 0;
    unordered_map<string, shared_ptr<DFA>> cache;
    double stats[8] = {};
    size_t cached_nodes = 0;
};
static thread_local string error;
static Part readpart(ifstream &f, const string &expect) {
    string name;
    int n, st, na;
    if (!(f >> name >> n >> st >> na) || name != expect || n <= 0 || st < 0 || st >= n || na < 0)
        throw runtime_error("invalid parts header: " + expect);
    Part p;
    p.start = st;
    p.fin.resize(n);
    p.arcs.resize(n);
    for (int &v : p.fin)
        if (!(f >> v) || v < 0 || v > INF)
            throw runtime_error("invalid final weight: " + expect);
    for (int i = 0; i < na; i++) {
        int s, a, b, w, t;
        if (!(f >> s >> a >> b >> w >> t) || s < 0 || s >= n || t < 0 || t >= n || a < 0 || a > 8 ||
            b < 0 || b > 8 || w < 0 || w >= INF)
            throw runtime_error("invalid primitive arc: " + expect);
        p.arcs[s][a].push_back({b, w, t});
    }
    return p;
}
static DFA linear(const string &s) {
    DFA d;
    d.len = s.size();
    d.start = 0;
    d.v.resize(s.size() + 1);
    for (int i = 0; i < (int)s.size(); i++) {
        int a = s[i] - '0';
        if (a < 0 || a > 8)
            throw runtime_error("invalid digit");
        d.v[i].to[a] = i + 1;
    }
    d.v.back().fin = 0;
    return d;
}
static DFA minimize(DFA d, const vector<vector<int>> &layers) {
    DFA o;
    o.len = d.len;
    vector<int> remap(d.v.size(), -1), pot(d.v.size(), INF);
    unordered_map<vector<int>, int, VH> memo;
    for (int lev = (int)layers.size() - 1; lev >= 0; lev--) {
        memo.clear();
        for (int s : layers[lev]) {
            const Node &n = d.v[s];
            int z = n.fin;
            for (int a = 0; a < 9; a++)
                if (n.to[a] >= 0 && pot[n.to[a]] < INF)
                    z = min(z, n.w[a] + pot[n.to[a]]);
            if (z >= INF)
                continue;
            pot[s] = z;
            vector<int> key(19, -1);
            key[0] = n.fin >= INF ? INF : n.fin - z;
            for (int a = 0; a < 9; a++)
                if (n.to[a] >= 0 && remap[n.to[a]] >= 0) {
                    key[1 + a * 2] = remap[n.to[a]];
                    key[2 + a * 2] = n.w[a] + pot[n.to[a]] - z;
                }
            auto it = memo.find(key);
            if (it != memo.end()) {
                remap[s] = it->second;
                continue;
            }
            int t = o.v.size();
            Node q;
            q.fin = key[0];
            for (int a = 0; a < 9; a++) {
                q.to[a] = key[1 + a * 2];
                q.w[a] = q.to[a] >= 0 ? key[2 + a * 2] : 0;
            }
            o.v.push_back(q);
            memo.emplace(move(key), t);
            remap[s] = t;
        }
    }
    if (d.start >= 0 && remap[d.start] >= 0) {
        o.start = remap[d.start];
        o.offset = d.offset + pot[d.start];
    }
    return o;
}
static DFA pull(Model *m, const Part &p, const DFA &a) {
    if (a.start < 0)
        return a;
    DFA d;
    d.offset = a.offset;
    d.len = a.len;
    d.start = 0;
    d.v.emplace_back();
    vector<vector<int>> layers(a.len + 1);
    layers[0].push_back(0);
    vector<vector<int>> cur = {{p.start, a.start, 0}};
    for (int lev = 0; lev <= a.len; lev++) {
        vector<vector<int>> nxt;
        unordered_map<vector<int>, int, VH> memo;
        for (size_t j = 0; j < cur.size(); j++) {
            const auto &subset = cur[j];
            int s = layers[lev][j];
            if (lev == a.len) {
                int z = INF;
                for (size_t k = 0; k < subset.size(); k += 3) {
                    int pf = p.fin[subset[k]], af = a.v[subset[k + 1]].fin;
                    if (pf < INF && af < INF)
                        z = min(z, subset[k + 2] + pf + af);
                }
                d.v[s].fin = z;
                continue;
            }
            for (int c = 0; c < 9; c++) {
                vector<pair<uint64_t, int>> cand;
                for (size_t k = 0; k < subset.size(); k += 3) {
                    int ps = subset[k], as = subset[k + 1], w = subset[k + 2];
                    const Node &an = a.v[as];
                    for (const Arc &r : p.arcs[ps][c])
                        if (an.to[r.out] >= 0) {
                            uint64_t key =
                                (uint64_t(uint32_t(r.to)) << 32) | uint32_t(an.to[r.out]);
                            cand.emplace_back(key, w + r.w + an.w[r.out]);
                        }
                }
                if (cand.empty())
                    continue;
                sort(cand.begin(), cand.end());
                vector<int> out;
                out.reserve(cand.size() * 3);
                int low = INF;
                for (size_t k = 0; k < cand.size();) {
                    uint64_t key = cand[k].first;
                    int val = cand[k].second;
                    size_t end = k + 1;
                    while (end < cand.size() && cand[end].first == key) {
                        val = min(val, cand[end].second);
                        end++;
                    }
                    out.push_back(int(key >> 32));
                    out.push_back(int(key & 0xffffffffu));
                    out.push_back(val);
                    low = min(low, val);
                    k = end;
                }
                for (size_t k = 2; k < out.size(); k += 3)
                    out[k] -= low;
                auto it = memo.find(out);
                int t;
                if (it != memo.end())
                    t = it->second;
                else {
                    t = d.v.size();
                    d.v.emplace_back();
                    layers[lev + 1].push_back(t);
                    memo.emplace(out, t);
                    nxt.push_back(move(out));
                }
                d.v[s].to[c] = t;
                d.v[s].w[c] = low;
            }
        }
        cur = move(nxt);
    }
    m->stats[0]++;
    m->stats[1] = max(m->stats[1], double(d.v.size()));
    DFA o = minimize(move(d), layers);
    m->stats[2] = max(m->stats[2], double(o.v.size()));
    return o;
}
static DFA half(const DFA &a) {
    if (a.start < 0)
        return a;
    DFA d = a;
    for (auto &n : d.v) {
        Node q;
        q.fin = n.fin;
        for (int c = 0; c <= 8; c++) {
            q.to[c] = n.to[wgd_model::double_cn(c)];
            q.w[c] = n.w[wgd_model::double_cn(c)];
        }
        n = q;
    }
    return d;
}
static shared_ptr<DFA> get(Model *m, const string &s, int br) {
    string key = char('a' + br + 1) + s;
    auto it = m->cache.find(key);
    if (it != m->cache.end()) {
        m->stats[3]++;
        return it->second;
    }
    DFA d;
    if (br < 0)
        d = pull(m, m->now, linear(s));
    else {
        string qkey = char('q' + br) + s;
        auto jt = m->cache.find(qkey);
        shared_ptr<DFA> q;
        if (jt != m->cache.end())
            q = jt->second;
        else {
            DFA v;
            if (br == 0)
                v = pull(m, m->lg, linear(s));
            else if (m->kind == 0) {
                get(m, s, br - 1);
                v = pull(m, m->lg, half(*m->cache.at(char('q' + br - 1) + s)));
            } else {
                get(m, s, 0);
                v = *m->cache.at('q' + s);
                for (int j = 0; j < br; j++)
                    v = half(v);
                if (m->kind == 1)
                    v = pull(m, m->lg, v);
            }
            q = make_shared<DFA>(move(v));
            m->cached_nodes += q->v.size();
            m->cache[qkey] = q;
        }
        d = pull(m, m->loh, *q);
    }
    auto r = make_shared<DFA>(move(d));
    m->cached_nodes += r->v.size();
    m->cache[key] = r;
    return r;
}
static int intersect(Model *m, const DFA &a, const DFA &b) {
    m->stats[4]++;
    if (a.start < 0 || b.start < 0)
        return INF;
    unordered_map<uint64_t, int> cur, nxt;
    cur[(uint64_t(uint32_t(a.start)) << 32) | uint32_t(b.start)] = a.offset + b.offset;
    for (int i = 0; i < a.len; i++) {
        nxt.clear();
        for (auto &kv : cur) {
            const auto &an = a.v[kv.first >> 32];
            const auto &bn = b.v[kv.first & 0xffffffffu];
            for (int c = 0; c < 9; c++)
                if (an.to[c] >= 0 && bn.to[c] >= 0) {
                    uint64_t k = (uint64_t(uint32_t(an.to[c])) << 32) | uint32_t(bn.to[c]);
                    int val = kv.second + an.w[c] + bn.w[c];
                    auto it = nxt.find(k);
                    if (it == nxt.end())
                        nxt[k] = val;
                    else if (val < it->second)
                        it->second = val;
                }
        }
        cur.swap(nxt);
        if (cur.empty())
            return INF;
    }
    int best = INF;
    for (auto &kv : cur) {
        int af = a.v[kv.first >> 32].fin, bf = b.v[kv.first & 0xffffffffu].fin;
        if (af < INF && bf < INF)
            best = min(best, kv.second + af + bf);
    }
    return best;
}
static vector<string> split(const char *s) {
    vector<string> out(1);
    for (; *s; s++) {
        if (*s == 'X')
            out.emplace_back();
        else {
            if (*s < '0' || *s > '8')
                throw runtime_error("expected digits 0..8 or X");
            out.back().push_back(*s);
        }
    }
    return out;
}
extern "C" {
const char *ie_error() {
    return error.c_str();
}
void *ie_load(const char *file) {
    try {
        auto m = make_unique<Model>();
        ifstream f(file);
        m->loh = readpart(f, "LOH");
        m->lg = readpart(f, "LG");
        m->now = readpart(f, "NOWGD");
        return m.release();
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}
void ie_set_model(void *h, int kind) {
    auto *m = static_cast<Model *>(h);
    m->kind = kind;
    m->cache.clear();
    m->cached_nodes = 0;
}
void ie_free(void *h) {
    delete static_cast<Model *>(h);
}
void ie_clear(void *h) {
    auto *m = static_cast<Model *>(h);
    m->cache.clear();
    m->cached_nodes = 0;
}
void ie_stats(void *h, double *out) {
    auto *m = static_cast<Model *>(h);
    for (int i = 0; i < 8; i++)
        out[i] = m->stats[i];
    out[6] = m->cached_nodes;
    out[7] = m->cache.size();
}
double ie_score(void *h, const char *x, const char *y, int directed) {
    try {
        auto *m = static_cast<Model *>(h);
        auto xs = split(x), ys = split(y);
        if (xs.size() != ys.size())
            return INFINITY;
        for (size_t i = 0; i < xs.size(); i++)
            if (xs[i].size() != ys[i].size())
                return INFINITY;
        int best = INF;
        vector<array<int, 25>> suffix(xs.size() +
                                      1); // value initialization sets every bound to zero
        for (int i = int(xs.size()) - 1; i >= 0; i--) {
            array<int, 16> lb{};
#ifdef IE_VARIATION_BOUND
            lb = factor_bound::chromosome(xs[i], ys[i]);
#endif
            array<int, 25> local{};
#ifdef IE_LOCAL_BOUND
            local = local_factor_bound::chromosome(xs[i], ys[i]);
#endif
            for (int a = -1; a <= IE_MAX_WGDS; a++)
                for (int b = -1; b <= IE_MAX_WGDS; b++) {
                    int j = 5 * (a + 1) + b + 1;
                    suffix[i][j] = suffix[i + 1][j] + max(local[j], lb[4 * max(a, 0) + max(b, 0)]);
                }
        }
        vector<array<int, 3>> branches;
        for (int a = -1; a <= IE_MAX_WGDS; a++)
            for (int b = -1; b <= IE_MAX_WGDS; b++)
                if (!directed || a == -1) {
                    int k = 5 * (a + 1) + b + 1;
                    branches.push_back({max(a, 0) + max(b, 0) + suffix[0][k], a, b});
                }
        sort(branches.begin(), branches.end());
        for (auto br : branches) {
            int a = br[1], b = br[2], k = 5 * (a + 1) + b + 1, val = max(a, 0) + max(b, 0);
            if (br[0] >= best)
                continue;
            m->stats[5]++;
            for (size_t i = 0; i < xs.size(); i++) {
                int v;
                if (directed) {
                    DFA lin = linear(xs[i]);
                    v = intersect(m, lin, *get(m, ys[i], b));
                } else if (a == b && a <= 0 && xs[i] == ys[i])
                    v = 0;
                else
                    v = intersect(m, *get(m, xs[i], a), *get(m, ys[i], b));
                if (v >= INF) {
                    val = INF;
                    break;
                }
                val += v;
                if (val + suffix[i + 1][k] >= best) {
                    val = INF;
                    break;
                }
            }
            best = min(best, val);
        }
        return best >= INF ? INFINITY : double(best);
    } catch (exception &e) {
        error = e.what();
        return -1;
    }
}
}
