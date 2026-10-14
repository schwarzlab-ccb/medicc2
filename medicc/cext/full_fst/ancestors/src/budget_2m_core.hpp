// Shared storage and primitive/DFA operations. The historical filename is NOT
// a two-million-state cap: the submitted/review configuration has no such cap.
// Read Node, DFA, readpart and minimize first. pull is older weighted-subset
// support, used by support paths rather than the direct upward/downward searches.
// Local adaptation from ../integer_engine/engine.cpp; includes X arcs.
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
using namespace std;

// Arc: output symbol, nonnegative integer weight, destination primitive state.
// Part is a small primitive transducer; arcs[state][input] lists possible arcs.
// Node is a deterministic acceptor state; to[c]==-1 means no transition.
// DFA owns all Nodes. offset stores the shared cost removed by weight pushing.
// vector and array clean up their owned buffers automatically on destruction.
struct Arc {
    int out, w, to;
};
struct Part {
    int start;
    vector<int> fin;
    vector<array<vector<Arc>, 10>> arcs;
};
#include "../../common/src/automata_types.hpp"
#include "../../common/src/openfst_language.hpp"
struct VH {
    size_t operator()(const vector<int> &v) const {
        size_t h = 0x9e3779b97f4a7c15ULL;
        for (int x : v)
            h ^= size_t(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};
// This support Model is different from AM (the ancestral model). shared_ptr
// denotes shared ownership of a DFA cache entry: the last owner deletes it.
// unordered_map is a hash table; its template arguments specify key/value types.
struct Model {
    Part loh, lg, now;
    int kind = 0, budget = INF;
    unordered_map<string, shared_ptr<DFA>> cache;
    double stats[8] = {};
    size_t cached_nodes = 0;
};
// Error text is per-thread. am_error returns a borrowed c_str() pointer;
// it stays valid until the next change to this thread's error string.
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
        if (!(f >> s >> a >> b >> w >> t) || s < 0 || s >= n || t < 0 || t >= n || a < 0 || a > 9 ||
            b < 0 || b > 9 || w < 0 || w >= INF)
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
        int a = s[i] == 'X' ? 9 : s[i] - '0';
        if (a < 0 || a > 9 || (s[i] != 'X' && s[i] > '8'))
            throw runtime_error("invalid digit");
        d.v[i].to[a] = i + 1;
    }
    d.v.back().fin = 0;
    return d;
}
// Process the layered DAG backward. Subtract each state's minimum suffix cost
// (weight pushing), then merge states with identical final/arc signatures.
// const Node& borrows one Node; unlike a value declaration, it does not copy it.
// This function owns d by value and constructs an independently owned result.
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
            for (int a = 0; a < 10; a++)
                if (n.to[a] >= 0 && pot[n.to[a]] < INF)
                    z = min(z, n.w[a] + pot[n.to[a]]);
            if (z >= INF)
                continue;
            pot[s] = z;
            vector<int> key(21, -1);
            key[0] = n.fin >= INF ? INF : n.fin - z;
            for (int a = 0; a < 10; a++)
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
            for (int a = 0; a < 10; a++) {
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
static vector<int> pull_potential(const Part &p, const DFA &a) {
    int np = p.fin.size();
    vector<int> pot(a.v.size() * np, INF);
    vector<vector<int>> ls(a.len + 1);
    vector<int> seen(a.v.size(), 0);
    ls[0] = {a.start};
    seen[a.start] = 1;
    for (int lev = 0; lev < a.len; lev++)
        for (int s : ls[lev])
            for (int c = 0; c < 10; c++) {
                int t = a.v[s].to[c];
                if (t >= 0 && !seen[t]) {
                    seen[t] = 1;
                    ls[lev + 1].push_back(t);
                }
            }
    for (int lev = a.len; lev >= 0; lev--)
        for (int s : ls[lev])
            for (int ps = 0; ps < np; ps++) {
                int z = p.fin[ps] + a.v[s].fin;
                for (int c = 0; c < 10; c++)
                    for (auto &r : p.arcs[ps][c]) {
                        int t = a.v[s].to[r.out];
                        if (t >= 0)
                            z = min(z, r.w + a.v[s].w[r.out] + pot[size_t(t) * np + r.to]);
                    }
                pot[size_t(s) * np + ps] = min(z, INF);
            }
    return pot;
}
static DFA pull(Model *m, const Part &p, const DFA &a) {
    if (a.start < 0)
        return a;
    vector<int> pot;
    if (m->budget < INF)
        pot = pull_potential(p, a);
    int np = p.fin.size();
    if (m->budget < INF && a.offset + pot[size_t(a.start) * np + p.start] > m->budget) {
        DFA e;
        e.len = a.len;
        return e;
    }
    vector<int> prefix(1, a.offset);
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
                        if (prefix[s] + subset[k + 2] + pf + af <= m->budget)
                            z = min(z, subset[k + 2] + pf + af);
                }
                d.v[s].fin = z;
                continue;
            }
            for (int c = 0; c < 10; c++) {
                vector<pair<uint64_t, int>> cand;
                for (size_t k = 0; k < subset.size(); k += 3) {
                    int ps = subset[k], as = subset[k + 1], w = subset[k + 2];
                    const Node &an = a.v[as];
                    for (const Arc &r : p.arcs[ps][c])
                        if (an.to[r.out] >= 0) {
                            uint64_t key =
                                (uint64_t(uint32_t(r.to)) << 32) | uint32_t(an.to[r.out]);
                            int value = w + r.w + an.w[r.out];
                            int h = pot.empty() ? 0 : pot[size_t(an.to[r.out]) * np + r.to];
                            if (prefix[s] + value + h <= m->budget)
                                cand.emplace_back(key, value);
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
                if (prefix.size() < d.v.size())
                    prefix.push_back(INF);
                prefix[t] = min(prefix[t], prefix[s] + low);
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
