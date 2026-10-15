// REVIEW START: bottom-up sibling join, with exact local prune_weight=0.
// Search primitive tuples first. After finding the best cost, retain EVERY
// optimal ancestor word. Only then determinize/minimize the surviving language.
// This preserves local OpenFST prune(0) semantics; it is not a claim of a global
// unpruned optimum for all internal profiles on a tree.
// Direct primitive tuples. Determinize only the globally optimal ancestor NFA.
#include "lazy_trace.hpp"
#include "../../common/src/feasible_path.hpp"
#include <chrono>
#include <sstream>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
static thread_local string join_trace;
using SearchClock = std::chrono::steady_clock;
static double elapsed(SearchClock::time_point t) { return std::chrono::duration<double>(SearchClock::now()-t).count(); }
// Shared across the outer tree tasks and inner WGD searches of ONE reconstruction.
// Waiting helper threads hold no permit. All expensive native work holds a permit.
struct SearchBudget {
    int capacity, used = 0, peak = 0;
    std::mutex mutex;
    std::condition_variable changed;
    explicit SearchBudget(int n) : capacity(n) {}
    void acquire() {
        std::unique_lock<std::mutex> lock(mutex);
        changed.wait(lock, [&] { return used < capacity; });
        peak = std::max(peak, ++used);
    }
    void release() {
        { std::lock_guard<std::mutex> lock(mutex); --used; }
        changed.notify_all();
    }
    template<class Remaining> bool acquire_helper(Remaining remaining) {
        std::unique_lock<std::mutex> lock(mutex);
        changed.wait(lock, [&] { return used < capacity || !remaining(); });
        if (!remaining()) return false;
        peak = std::max(peak, ++used);
        return true;
    }
};
static thread_local SearchBudget *worker_budget = nullptr;
#include "generated/chromosome_bound.hpp"
#include "../../common/src/boost_search.hpp"
namespace direct_nfa {
struct Edge {
    uint64_t to;
    int w, h, out;
};
// One possible history branch for one child. child and parts are BORROWED
// pointers: Side never deletes them. The caller keeps the DFA and AM model alive.
// vector<const Part*> owns its pointer array, but NOT the pointed-to Parts.
// const applies to the Part: this search may read, but cannot change, its arcs.
struct Side {
    const DFA *child;
    vector<const Part *> parts;
    vector<uint32_t> radix;
    int wg;
    uint32_t startpart = 0;
    shared_ptr<const vector<int>> pot, chrom_pot;
    int np;
    Side(const DFA &d, vector<const Part *> chain, int k, int multiple,
         shared_ptr<const vector<int>> shared_pot = nullptr,
         shared_ptr<const vector<int>> shared_chrom = nullptr, bool with_chrom = true)
        : child(&d), parts(move(chain)), wg(k) {
        // Mixed-radix packing: q0 + n0*q1 + n0*n1*q2 + ... stores primitive
        // states in one integer. radix[i] is the multiplier for primitive i.
        // This representation assumes the supported small primitives fit in 32 bits;
        // changing the alphabet/model requires revisiting that representation.
        if (parts.empty() || parts.size() > 8) throw runtime_error("unsupported primitive depth");
        uint32_t mul = 1;
        for (auto *p : parts) {
            if (p->fin.empty() || p->fin.size() > UINT32_MAX / mul)
                throw runtime_error("primitive tuple exceeds 32-bit packing");
            radix.push_back(mul);
            startpart += p->start * mul;
            mul *= p->fin.size();
        }
        // pot[child_state * np + last_primitive_state] is an exact relaxed suffix
        // cost. It ignores all earlier primitives, so it is an admissible lower bound.
        // Final weights matter too; chromosome X is handled as an ordinary explicit arc.
        np = parts.back()->fin.size();
        if (shared_pot) { pot = move(shared_pot); chrom_pot = move(shared_chrom); return; }
        auto values = make_shared<vector<int>>(d.v.size() * np, INF);
        auto partial = with_chrom ? make_shared<vector<int>>(d.v.size() * np, INF) : nullptr;
        auto ls = levels(d);
        auto &p = *parts.back();
        for (int lev = d.len; lev >= 0; lev--)
            for (int s : ls[lev])
                for (int q = 0; q < np; q++) {
                    int z = min(INF, p.fin[q] + d.v[s].fin), cz = z;
                    for (int c = 0; c < 10; c++)
                        if (c == 9 || c % multiple == 0)
                            for (auto &e : p.arcs[q][c]) {
                                int t = d.v[s].to[e.out];
                                if (t >= 0) {
                                    z = min(z, e.w + d.v[s].w[e.out] + (*values)[size_t(t) * np + e.to]);
                                    if (partial) cz = min(cz, e.w + d.v[s].w[e.out] +
                                        (c == 9 ? 0 : (*partial)[size_t(t) * np + e.to]));
                                }
                            }
                    (*values)[size_t(s) * np + q] = min(z, INF);
                    if (partial) (*partial)[size_t(s) * np + q] = min(cz, INF);
                }
        pot = move(values);
        chrom_pot = move(partial);
    }
    // State key: high 32 bits = primitive tuple, low 32 bits = child DFA state.
    // << shifts bits; | joins the disjoint bit fields; uint32_t extracts low bits.
    // These are numeric state IDs, NOT memory addresses.
    uint64_t start() const {
        return uint64_t(startpart) << 32 | uint32_t(child->start);
    }
    int h(uint64_t s) const {
        uint32_t p = s >> 32;
        int q = (p / radix.back()) % np;
        return (*pot)[size_t(uint32_t(s)) * np + q];
    }
    int chrom_h(uint64_t s) const {
        uint32_t p = s >> 32;
        int q = (p / radix.back()) % np;
        return (*chrom_pot)[size_t(uint32_t(s)) * np + q];
    }
    int final(uint64_t s) const {
        uint32_t p = s >> 32;
        int val = child->v[uint32_t(s)].fin;
        for (size_t j = 0; j < parts.size(); j++)
            val = min(INF, val + parts[j]->fin[(p / radix[j]) % parts[j]->fin.size()]);
        return val;
    }
    // Stack-owned lookahead tables for one column. array<T,N> is a fixed-size
    // C++ array; vector<T> is a resizable contiguous array with automatic cleanup.
    struct Prepared {
        int ds;
        array<int, 8> ps{};
        array<array<int, 10>, 9> lower;
    };
    Prepared prepare(uint64_t state) const {
        Prepared z{};
        uint32_t code = state >> 32;
        int ds = uint32_t(state), n = parts.size();
        z.ds = ds;
        auto &ps = z.ps;
        auto &lower = z.lower;
        for (int i = 0; i < n; i++)
            ps[i] = (code / radix[i]) % parts[i]->fin.size();
        // Per-column lookahead. Retain the last primitive's exact suffix relaxation.
        lower[n].fill(INF);
        for (int i = n - 1; i >= 0; i--)
            for (int c = 0; c < 10; c++) {
                int val = INF;
                for (auto &a : parts[i]->arcs[ps[i]][c]) {
                    if (i == n - 1) {
                        int t = child->v[ds].to[a.out];
                        if (t >= 0)
                            val =
                                min(val, a.w + child->v[ds].w[a.out] + (*pot)[size_t(t) * np + a.to]);
                    } else
                        val = min(val, a.w + lower[i + 1][a.out]);
                }
                lower[i][c] = min(val, INF);
            }
        return z;
    }
    vector<Edge> edges(const Prepared &z, int input, int cap) const {
        vector<Edge> out;
        if (cap < 0)
            return out;
        int ds = z.ds, n = parts.size();
        const auto &ps = z.ps;
        const auto &lower = z.lower;
        // Recursive generic lambda. auto&& accepts the lambda itself as a parameter
        // so it can call visit(visit, ...). All captures are borrowed and remain local.
        // Enumerate compatible outputs through the small primitives for ONE column;
        // do not build a composed transducer for all possible strings.
        auto visit = [&](auto &&visit, int i, int label, uint32_t next, int w) -> void {
            if (w + lower[i][label] > cap)
                return;
            for (auto &a : parts[i]->arcs[ps[i]][label]) {
                int nw = w + a.w;
                uint32_t nc = next + radix[i] * a.to;
                if (i + 1 == n) {
                    int t = child->v[ds].to[a.out];
                    if (t < 0)
                        continue;
                    int weight = nw + child->v[ds].w[a.out];
                    int tail = (*pot)[size_t(t) * np + a.to];
                    if (weight + tail <= cap)
                        out.push_back({uint64_t(nc) << 32 | uint32_t(t), weight, tail, a.out});
                } else
                    visit(visit, i + 1, a.out, nc, nw);
            }
        };
        if (lower[0][input] <= cap)
            visit(visit, 0, input, 0, 0);
        // [] means the comparison lambda captures nothing. tie(a.to,a.w) builds
        // a tuple of references for lexicographic sorting. For identical next states,
        // only the cheapest path is needed. The upward output label is the common
        // input symbol; downward traceback needs only one optimal emitted symbol.
        sort(out.begin(), out.end(),
             [](const Edge &a, const Edge &b) { return tie(a.to, a.w) < tie(b.to, b.w); });
        size_t w = 0;
        for (size_t k = 0; k < out.size();) {
            size_t e = k + 1;
            while (e < out.size() && out[e].to == out[k].to)
                e++;
            out[w++] = out[k];
            k = e;
        }
        out.resize(w);
        sort(out.begin(), out.end(),
             [](const Edge &a, const Edge &b) { return a.w + a.h < b.w + b.h; });
        return out;
    }
};
struct StateKey {
    uint64_t a, b;
    bool operator==(const StateKey &o) const {
        return a == o.a && b == o.b;
    }
};
struct StateKeyHash {
    size_t operator()(const StateKey &k) const {
        size_t x = k.a ^ (k.b + 0x9e3779b97f4a7c15ULL + (k.a << 6) + (k.a >> 2));
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        return x;
    }
};
struct SearchState {
    StateKey key;
    int level, g = INF, tail = INF;
    bool expanded = false;
};
// Open-addressed state index. Keys remain owned by nodes; the table stores
// only integer IDs, so vector reallocations cannot invalidate its entries.
// At <=70% occupancy every unsuccessful lookup reaches an empty slot.
class StateIndex {
    const vector<SearchState> &nodes;
    vector<int> slots;
    size_t used = 0;
    void place(int id) {
        size_t pos = StateKeyHash{}(nodes[id].key) & (slots.size() - 1);
        while (slots[pos] >= 0) pos = (pos + 1) & (slots.size() - 1);
        slots[pos] = id;
    }
public:
    explicit StateIndex(const vector<SearchState> &n) : nodes(n), slots(256, -1) {}
    int find(const StateKey &key) const {
        size_t pos = StateKeyHash{}(key) & (slots.size() - 1);
        while (slots[pos] >= 0) {
            int id = slots[pos];
            if (nodes[id].key == key) return id;
            pos = (pos + 1) & (slots.size() - 1);
        }
        return -1;
    }
    void add(int id) {
        if ((used + 1) * 10 > slots.size() * 7) {
            vector<int> old;
            old.swap(slots);
            slots.assign(old.size() * 2, -1);
            for (int previous : old) if (previous >= 0) place(previous);
        }
        place(id);
        ++used;
    }
};
// RAII: the destructor ~Stats runs automatically at scope exit, even on a
// C++ exception. v is a borrowed output array; the destructor updates counters
// but must not delete that array.
struct Stats {
    double *v;
    size_t expanded = 0, states = 0, generated = 0, tight = 0, det = 0;
    size_t relaxation_skips = 0;
    double forward_seconds = 0, backward_seconds = 0, output_seconds = 0;
    ~Stats() {
        v[0] += expanded;
        v[1] = max(v[1], double(states));
        v[2] += generated;
        v[3] = max(v[3], double(tight));
        v[4] = max(v[4], double(det));
    }
};
// The tight NFA now contains only paths of one common optimal cost. Therefore
// subsets are sets of state IDs, not weighted residual-cost subsets. Store that
// common cost in d.offset and keep the remaining transition weights zero.
// move(...) transfers a vector/DFA buffer to its next owner instead of copying.
static DFA determinize(const vector<Tight> &n, int start, int len, int offset, Stats &stats,
                       int limit) {
    size_t states=0;
    auto d=library_language::determinized(n,start,len,offset,states);
    stats.det=max(stats.det,states);
    auto layers=levels(d);
    return minimize(move(d),layers);
}
// Joint state = left Side state plus right Side state. Both consume the same
// ancestor symbol. upper is a cost of a FEASIBLE complete witness, not a beam.
// Only f > cap can be discarded here: equality may encode another optimal word.
static DFA search(Side &a, Side &b, const vector<int> &mask, int upper, int limit, Stats &stats, const vector<int> &future,
                  std::atomic<int> *incumbent = nullptr) {
    DFA empty;
    empty.len = a.child->len;
    int offset = a.child->offset + b.child->offset + a.wg + b.wg;
    int cap = upper >= INF ? INF : upper - offset;
    int h0 = max(a.h(a.start()) + b.h(b.start()),
                 a.chrom_h(a.start()) + b.chrom_h(b.start()) + future[0]);
    if (h0 >= INF || h0 > cap)
        return empty;
    vector<SearchState> nodes;
    vector<int> heuristic;
    StateIndex ids(nodes);
    auto insert = [&](StateKey k, int lev) {
        int found = ids.find(k);
        if (found >= 0) return found;
        int s = nodes.size();
        nodes.push_back({k, lev});
        heuristic.push_back(INF);
        ids.add(s);
        stats.states = max(stats.states, nodes.size());
        return s;
    };
    int root = insert({a.start(), b.start()}, 0);
    nodes[root].g = 0;
    heuristic[root] = h0;
    bool found = false;
    int found_cost = INF;
    // IDs in nodes remain stable if its vector reallocates. Capture key/g/level
    // by value before inserting new states. Child Side objects stay alive throughout.
    // No product arcs or intermediate subsets are cached. Regenerate from tuples.
    auto transitions = [&](StateKey key, int lev, int remaining, auto &&fn) {
        auto ap = a.prepare(key.a), bp = b.prepare(key.b);
        for (int c = 0; c < 10; c++)
            if (mask[lev] & (1 << c)) {
                // Both sides must emit this SAME ancestor symbol. Their
                // one-column lookaheads bound disjoint branch costs.
                if (ap.lower[0][c] + bp.lower[0][c] > remaining) continue;
                auto left = a.edges(ap, c, remaining - bp.lower[0][c]);
                if (left.empty())
                    continue;
                auto right = b.edges(bp, c, remaining - left.front().w - left.front().h);
                if (right.empty())
                    continue;
                for (auto &x : left) {
                    int room = remaining - x.w - x.h;
                    for (auto &y : right) {
                        if (y.w + y.h > room)
                            break;
                        // Current-chromosome costs and complete later-chromosome
                        // costs are disjoint. Never add this to the full old bound.
                        int tail = max(x.h + y.h,
                            a.chrom_h(x.to) + b.chrom_h(y.to) + future[lev + 1]);
                        if (x.w + y.w + tail > remaining) continue;
                        stats.generated++;
                        fn(c, StateKey{x.to, y.to}, x.w + y.w, tail);
                    }
                }
            }
    };
    auto forward_start = SearchClock::now();
    library_search::run(root, nodes, [&](int s, vector<library_search::Edge> &edges) {
        int g = nodes[s].g;
        if (incumbent) cap = min(cap, incumbent->load(std::memory_order_relaxed) - offset);
        if (g + heuristic[s] > cap) throw library_search::Stop{};
        stats.expanded++;
        nodes[s].expanded = true;
        StateKey key = nodes[s].key;
        int lev = nodes[s].level;
        if (lev == empty.len) {
            int fin = a.final(key.a) + b.final(key.b);
            if (fin < INF && g + fin <= cap) {
                cap = g + fin;
                found = true;
                found_cost = cap;
                if (incumbent) {
                    int total = offset + cap, previous = incumbent->load(std::memory_order_relaxed);
                    while (total < previous && !incumbent->compare_exchange_weak(
                               previous, total, std::memory_order_relaxed)) {}
                }
            }
            return;
        }
        transitions(key, lev, cap - g, [&](int c, StateKey k, int w, int tail) {
            int t = insert(k, lev + 1);
            heuristic[t] = tail;
            // Boost would reject this same non-improving distance. Avoid
            // materializing the edge; the backward pass regenerates ALL labels
            // and ties independently from the retained states.
            if (g + w < nodes[t].g) edges.push_back({s,t,w});
            else ++stats.relaxation_skips;
        });
    }, [&](int s) {return library_search::Rank{heuristic[s],0,s};});
    stats.forward_seconds += elapsed(forward_start);
    auto backward_start = SearchClock::now();
    // Backward pass regenerates transitions instead of retaining all explored arcs.
    // tail is exact on the explored region. All potentially optimal states survived
    // the admissible A* cutoff, including equal-cost states.
    if (!found || found_cost > cap)
        return empty;
    vector<vector<int>> layers(empty.len + 1);
    for (int s = 0; s < int(nodes.size()); s++)
        if (nodes[s].expanded && nodes[s].g + a.h(nodes[s].key.a)
                + b.h(nodes[s].key.b) <= cap)
            layers[nodes[s].level].push_back(s);
    for (int lev = empty.len; lev >= 0; lev--)
        for (int s : layers[lev]) {
            auto key = nodes[s].key;
            int g = nodes[s].g;
            if (lev == empty.len) {
                nodes[s].tail = min(INF, a.final(key.a) + b.final(key.b));
                continue;
            }
            transitions(key, lev, cap - g, [&](int, StateKey k, int w, int) {
                int target = ids.find(k);
                if (target >= 0)
                    nodes[s].tail = min(nodes[s].tail, w + nodes[target].tail);
            });
        }
    if (nodes[root].tail != cap)
        throw runtime_error("direct NFA backward optimum mismatch");
    // A state is tight iff g + tail == optimum. An arc is tight iff
    // g(source) + arc_cost + tail(target) == optimum. Retain the union of ALL such
    // paths. This is why the next tree join still sees all locally optimal ancestors.
    vector<int> remap(nodes.size(), -1);
    vector<Tight> kept;
    for (int s = 0; s < int(nodes.size()); s++)
        if (nodes[s].g + nodes[s].tail == cap) {
            remap[s] = kept.size();
            kept.emplace_back();
        }
    stats.tight = max(stats.tight, kept.size());
    for (int s = 0; s < int(nodes.size()); s++)
        if (remap[s] >= 0) {
            int lev = nodes[s].level;
            auto key = nodes[s].key;
            int g = nodes[s].g;
            int dest = remap[s];
            if (lev == empty.len) {
                kept[dest].final = true;
                continue;
            }
            transitions(key, lev, cap - g, [&](int c, StateKey k, int w, int) {
                int t = ids.find(k);
                if (t >= 0) {
                    if (g + w + nodes[t].tail == cap) {
                        if (remap[t] < 0)
                            throw runtime_error("missing tight NFA target");
                        kept[dest].to[c].push_back(remap[t]);
                    }
                }
            });
            for (auto &v : kept[dest].to) {
                sort(v.begin(), v.end());
                v.erase(unique(v.begin(), v.end()), v.end());
            }
        }
    stats.backward_seconds += elapsed(backward_start);
    auto output_start = SearchClock::now();
    auto result = determinize(kept, remap[root], empty.len, offset + cap, stats, limit);
    stats.output_seconds += elapsed(output_start);
    return result;
}
// unique_ptr<Side> owns one heap-allocated Side. make_unique constructs it;
// the Side is deleted automatically when its unique_ptr leaves the vector/scope.
// Its borrowed Part/DFA pointers are still NOT owned and are not deleted.
// Full L_LOH cascade: A (D A)^k, A=L_LOH*G, k=0..2.
// The certified A*A=A identity collapses skipped optional W1 slots.
// Mandatory D has zero arc cost; charge k once for the entire genome.
static vector<unique_ptr<Side>> branches(AM *m, const DFA &d, int kind, bool with_chrom = true) {
    if (kind != 0) throw runtime_error("L_LOH pilot supports only the full cascade");
    vector<unique_ptr<Side>> out;
    for (int k = 0; k <= ANCESTRAL_MAX_WGDS; k++) {
        vector<const Part *> chain = {&m->a};
        for (int j = 0; j < k; j++) {
            chain.push_back(&m->d);
            chain.push_back(&m->a);
        }
        // Fused A receives even CN directly from the last mandatory D.
        // Its internal loss/gain histories remain unrestricted.
        out.push_back(make_unique<Side>(d, chain, k, k > 0 ? 2 : 1,
            k == 2 ? out[1]->pot : nullptr,
            k == 2 ? out[1]->chrom_pot : nullptr, with_chrom));
    }
    return out;
}

} // namespace direct_nfa
// Python-facing join: h/l/r are borrowed opaque pointers (model/left/right).
// Returns a NEW heap DFA. Ownership passes to Message; am_drop must delete it once.
// stat points to five doubles owned by Python. kind=0 is full LG; limit=0 is uncapped.
// nullptr means an exception; am_error contains a thread-local diagnostic.
extern "C" void *am_nfa_join(void *h, void *l, void *r, int upper, int kind, int limit,
                             double *stat) {
    try {
        using namespace direct_nfa;
        auto *m = (AM *)h;
        auto &a = *(DFA *)l;
        auto &b = *(DFA *)r;
        auto setup_start = SearchClock::now();
        join_trace.clear();
        if (kind < 0 || kind > 1 || limit != 0)
            throw runtime_error("invalid NFA configuration");
        if (a.len != b.len)
            throw runtime_error("length mismatch");
        if (a.start < 0 || b.start < 0) {
            auto *d = new DFA();
            d->len = a.len;
            return d;
        }
        auto left = branches(m, a, kind), right = branches(m, b, kind);
        auto lm = lazy_cascade::layout(a), rm = lazy_cascade::layout(b);
        vector<int> mask(a.len);
        for (int i = 0; i < a.len; i++) {
            if ((lm[i] & 511) && (rm[i] & 511))
                mask[i] = 510 | ((lm[i] & 1) && (rm[i] & 1) ? 1 : 0);
            if ((lm[i] & 512) && (rm[i] & 512))
                mask[i] |= 512;
        }
        vector<array<int, 3>> order;
        for (int i = 0; i < int(left.size()); i++)
            for (int j = 0; j < int(right.size()); j++) {
                auto &x = *left[i];
                auto &y = *right[j];
                int lb = a.offset + b.offset + x.wg + y.wg + x.h(x.start()) + y.h(y.start());
                order.push_back({lb, i, j});
            }
        sort(order.begin(), order.end());
        // Keep all branches tying at the global local-join optimum; never discard a
        // branch just because its lower bound EQUALS the current best cost.
        double setup_seconds = elapsed(setup_start);
        std::atomic<int> best{upper}, next{0};
        std::atomic<bool> failed{false};
        vector<DFA> results(order.size());
        vector<array<double, 5>> counters(order.size());
        vector<string> traces(order.size());
        vector<std::exception_ptr> errors(order.size());
        auto *budget = worker_budget;
        auto run = [&](int index) {
            try {
                const auto item = order[index];
                int bound = best.load(std::memory_order_relaxed);
                if (item[0] > bound) return;
                Stats stats{counters[index].data()};
                auto branch_start = SearchClock::now();
                auto future = chromosome_bound::future(lm, rm, 3*item[1] + item[2]);
                // A complete walk through both child languages supplies a safe
                // incumbent. It cannot remove an optimal or equally optimal word.
                auto &x = *left[item[1]], &y = *right[item[2]];
                int offset = a.offset + b.offset + x.wg + y.wg;
                int initial_bound = bound;
                bound = offset + feasible_path::cost(x, y, mask, bound - offset);
                int previous = best.load(std::memory_order_relaxed);
                while (bound < previous && !best.compare_exchange_weak(previous, bound)) {}
                bound = min(bound, best.load(std::memory_order_relaxed));
                double witness_seconds = elapsed(branch_start);
                auto d = search(x, y, mask, bound, limit, stats, future, &best);
                std::ostringstream trace;
                trace << "{\"left_wgd\":" << item[1] << ",\"right_wgd\":" << item[2]
                      << ",\"seconds\":" << elapsed(branch_start)
                      << ",\"initial_upper\":" << initial_bound
                      << ",\"witness_seconds\":" << witness_seconds
                      << ",\"upper\":" << bound << ",\"lower\":" << item[0]
                      << ",\"score\":" << (d.start < 0 ? -1 : d.offset)
                      << ",\"expansions\":" << stats.expanded << ",\"states\":" << stats.states
                      << ",\"generated\":" << stats.generated
                      << ",\"relaxation_skips\":" << stats.relaxation_skips
                      << ",\"forward\":" << stats.forward_seconds
                      << ",\"backward\":" << stats.backward_seconds
                      << ",\"output\":" << stats.output_seconds << "}";
                traces[index] = trace.str();
                results[index] = move(d);
            } catch (...) {
                errors[index] = std::current_exception();
                failed.store(true);
            }
        };
        auto remaining = [&] { return !failed.load() && next.load() < int(order.size()); };
        auto helper = [&] {
            while (budget->acquire_helper(remaining)) {
                int index = next.fetch_add(1);
                if (index < int(order.size()) && !failed.load()) run(index);
                budget->release();
            }
        };
        vector<std::thread> helpers;
        if (budget) {
            helpers.reserve(budget->capacity - 1);
            for (int i = 1; i < budget->capacity; i++) {
                // If the OS refuses an additional thread, the owner still handles
                // remaining work. Already-created threads must always be joined.
                try { helpers.emplace_back(helper); } catch (const std::system_error &) { break; }
            }
        }
        while (remaining()) {
            int index = next.fetch_add(1);
            if (index < int(order.size())) run(index);
        }
        if (budget) {
            budget->changed.notify_all();
            // Do not occupy a CPU-work permit while waiting for helper completion.
            // Reacquire before returning to the Python worker's scoped permit.
            budget->release();
            for (auto &thread : helpers) thread.join();
            budget->acquire();
        }
        for (auto &err : errors) if (err) std::rethrow_exception(err);
        vector<DFA> kept;
        std::ostringstream trace;
        trace << "{\"setup_seconds\":" << setup_seconds << ",\"initial_upper\":" << upper
              << ",\"helper_threads\":" << helpers.size() << ",\"branches\":[";
        bool first_trace = true;
        for (size_t index = 0; index < results.size(); index++) {
            for (int i = 0; i < 5; i++) {
                if (i == 0 || i == 2) stat[i] += counters[index][i];
                else stat[i] = max(stat[i], counters[index][i]);
            }
            if (!traces[index].empty()) {
                if (!first_trace) trace << ",";
                first_trace = false;
                trace << traces[index];
            }
            auto &d = results[index];
            if (d.start >= 0 && d.offset == best.load()) kept.push_back(move(d));
        }
        trace << "]}";
        join_trace = trace.str();
        return new DFA(unite(kept));
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}

// Each OS thread has its own output buffer, so independent workers do not race.
// c_str() below is borrowed, not a malloc/new allocation to be freed by the caller.
static thread_local string selected_profile;
extern "C" const char *am_choose_profile(void *h, int *cost) {
    try {
        auto &d = *(DFA *)h;
        if (d.start < 0)
            throw runtime_error("No feasible ancestor");
        int s = d.start;
        selected_profile.clear();
        selected_profile.reserve(d.len);
        *cost = d.offset;
        for (int i = 0; i <= d.len; i++) {
            if (d.v[s].fin == 0)
                return selected_profile.c_str();
            bool found = false;
            for (int c = 0; c < 10; c++)
                if (d.v[s].to[c] >= 0 && d.v[s].w[c] == 0) {
                    selected_profile.push_back(c == 9 ? 'X' : char('0' + c));
                    s = d.v[s].to[c];
                    found = true;
                    break;
                }
            if (!found)
                throw runtime_error("Dead state in retained message");
        }
        throw runtime_error("Path exceeds message length");
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}

extern "C" int am_max_wgds() {
    return ANCESTRAL_MAX_WGDS;
}

extern "C" int am_lloh_model() { return 1; }

extern "C" int am_fused_model() { return 1; }

extern "C" const char *am_last_join_trace() { return join_trace.c_str(); }

extern "C" void *am_budget_new(int workers) {
    try {
        if (workers < 1) throw runtime_error("invalid shared worker budget");
        return new SearchBudget(workers);
    } catch (exception &e) { error = e.what(); return nullptr; }
}
extern "C" void am_budget_free(void *b) { delete static_cast<SearchBudget *>(b); }
extern "C" int am_worker_enter(void *b) {
    if (!b || worker_budget) { error = "invalid or nested worker permit"; return 0; }
    auto *budget = static_cast<SearchBudget *>(b);
    budget->acquire(); worker_budget = budget; return 1;
}
extern "C" void am_worker_leave() {
    auto *budget = worker_budget;
    worker_budget = nullptr;
    if (budget) budget->release();
}
extern "C" int am_budget_peak(void *b) {
    auto *budget = static_cast<SearchBudget *>(b);
    std::lock_guard<std::mutex> lock(budget->mutex);
    return budget->peak;
}
