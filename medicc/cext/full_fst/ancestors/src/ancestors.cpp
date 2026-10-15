// REVIEW START: top-down reconstruction of ONE profile.
// Given a fixed parent string p and a retained child message M, minimize
//                  T(p -> c) + M(c).
// The state stores primitive FST states and a child-DFA state, not an eager T.
// No OpenFST composition or weighted determinization happens in this search.
// See nfa_engine.cpp for the bottom-up operation, which MUST retain all ties.
// One-sided primitive tuple A*: choose one optimal retained child profile.
#include "nfa_engine.hpp"
namespace direct_trace {
// prev is an INTEGER INDEX into nodes, not a pointer. vector growth can move
// the storage; indices remain valid while pointers/references to its elements
// would not. label stores the emitted child symbol for traceback.
struct Node {
    uint64_t key;
    int level, g = INF, prev = -1, label = -1;
};
struct Result {
    int cost = INF;
    string word;
};
// Side& is a borrowed reference (like a non-null C pointer used with dot syntax).
// const string& borrows the parent without copying it and forbids modification.
// double* stats points to three caller-owned doubles. No ownership is transferred.
// Result is returned by value; its std::string manages its own allocated buffer.
static Result search(direct_nfa::Side &a, const string &parent, int upper, double *stats) {
    Result result;
    int offset = a.child->offset + a.wg;
    int cap = upper >= INF ? INF : upper - offset;
    // Independent first-primitive relaxation conditioned on the actual parent.
    // It and the last-primitive/child potential concern disjoint primitive costs.
    // parts.front() is a pointer borrowed from the model. * dereferences it;
    // p is an alias for that Part, not a copy. The model outlives this call.
    const Part &p = *a.parts.front();
    int np = p.fin.size(), len = parent.size();
    vector<int> first(size_t(len + 1) * np, INF);
    for (int q = 0; q < np; q++)
        first[size_t(len) * np + q] = p.fin[q];
    // The zero-WGD fused branch uses the same sole primitive in a.h.
    // Its separate prefix relaxation would be discarded, so avoid computing it.
    for (int i = a.parts.size() > 1 ? len - 1 : -1; i >= 0; i--) {
        int c = parent[i] == 'X' ? 9 : parent[i] - '0';
        for (int q = 0; q < np; q++)
            for (const auto &e : p.arcs[q][c])
                first[size_t(i) * np + q] =
                    min(first[size_t(i) * np + q], e.w + first[size_t(i + 1) * np + e.to]);
    }
    // [&] creates a lambda (local function object) borrowing surrounding variables.
    // The heuristic adds two disjoint primitive costs: first primitive conditioned
    // on the parent's suffix, and last primitive plus child-message suffix.
    // Ignoring the middle primitives/compatibility can only lower the true cost.
    auto h = [&](uint64_t key, int lev) {
        int q = uint32_t(key >> 32) % np;
        // A zero-WGD branch has only one fused primitive: avoid double counting.
        return min(INF, a.h(key) + (a.parts.size() > 1 ? first[size_t(lev) * np + q] : 0));
    };
    int h0 = h(a.start(), 0);
    if (h0 >= INF || h0 > cap)
        return result;
    vector<Node> nodes;
    vector<int> heuristic;
    unordered_map<uint64_t, int> ids;
    // unordered_map is a hash table: packed state -> integer index in nodes.
    // Child states imply the column because messages are fixed-length layered DAGs.
    // Never keep a Node reference across insert(): push_back may move the vector.
    auto insert = [&](uint64_t key, int lev) {
        auto it = ids.find(key);
        if (it != ids.end())
            return it->second;
        int s = nodes.size();
        nodes.push_back({key, lev});
        heuristic.push_back(INF);
        ids.emplace(key, s);
        stats[1] = max(stats[1], double(nodes.size()));
        return s;
    };
    int root = insert(a.start(), 0);
    nodes[root].g = 0;
    // At equal f prefer deeper prefixes: only one optimum path is required.
    // The Boost.Graph adapter preserves (estimated total, negative depth,
    // cost so far, node index) ordering and records predecessors on relaxation.
    heuristic[root] = h0;
    int winner = -1;
    library_search::run(root, nodes, [&](int s, vector<library_search::Edge> &edges) {
        int g = nodes[s].g, f = g + heuristic[s];
        if (f > cap || (winner >= 0 && f >= cap)) throw library_search::Stop{};
        stats[0]++;
        auto key = nodes[s].key;
        int lev = nodes[s].level;
        if (lev == len) {
            int fin = a.final(key);
            if (fin < INF && g + fin <= cap) {cap = g + fin; winner = s;}
            return;
        }
        int input = parent[lev] == 'X' ? 9 : parent[lev] - '0';
        auto prep = a.prepare(key);
        for (const auto &e : a.edges(prep, input, cap - g)) {
            stats[2]++;
            int tail = h(e.to, lev + 1);
            if (tail >= INF || g + e.w + tail > cap) continue;
            int t = insert(e.to, lev + 1);
            heuristic[t] = tail;
            edges.push_back({s,t,e.w,e.out});
        }
    }, [&](int s) {return library_search::Rank{heuristic[s],0,s,-nodes[s].level};},
       [&](library_search::Edge e) {nodes[e.to].prev=e.from;nodes[e.to].label=e.label;});
    // Every predecessor link refers to an actually explored feasible transition.
    // A numerical bound alone is never treated as a recovered path.
    if (winner < 0)
        return result;
    result.cost = offset + cap;
    result.word.resize(len);
    int s = winner;
    for (int i = len - 1; i >= 0; i--) {
        int c = nodes[s].label;
        if (c < 0)
            throw runtime_error("missing direct traceback predecessor");
        result.word[i] = c == 9 ? 'X' : char('0' + c);
        s = nodes[s].prev;
    }
    if (s != root)
        throw runtime_error("invalid direct traceback length");
    return result;
}
} // namespace direct_trace
// C ABI entry point used by Python ctypes. extern "C" prevents C++ name mangling.
// h: borrowed AM* model; c: borrowed DFA* child message; neither is freed here.
// parent: borrowed NUL-terminated input string, copied to x for this call.
// score: caller-owned output integer. stats: caller-owned array of three doubles.
// Return: BORROWED thread-local string buffer, valid until this thread's next
// resolving call. Python c_char_p copies it immediately. Do NOT free it in Python.
// Exceptions are caught here: nullptr signals an error, retrievable by am_error.
extern "C" const char *am_direct_resolve(void *h, const char *parent, void *c, int upper,
                                         int *score, double *stats) {
    try {
        // Cast opaque C handles back to their real types. auto* infers AM*;
        // auto& child aliases the existing DFA. This does not allocate either object.
        auto *m = (AM *)h;
        auto &child = *(DFA *)c;
        string x = parent;
        if (int(x.size()) != child.len)
            throw runtime_error("length mismatch");
        for (char ch : x)
            if (ch != 'X' && (ch < '0' || ch > '8'))
                throw runtime_error("invalid parent symbol");
        if (child.start < 0)
            throw runtime_error("No feasible conditional ancestor");
        auto branches = direct_nfa::branches(m, child, 0,
            /*with_chrom=*/false, /*compact_wgd=*/false);
        vector<pair<int, int>> order;
        for (int i = 0; i < int(branches.size()); i++) {
            auto &a = *branches[i];
            order.push_back({child.offset + a.wg + a.h(a.start()), i});
        }
        sort(order.begin(), order.end());
        int best = upper;
        bool found = false;
        string word;
        for (auto item : order) {
            if (item.first > best || (found && item.first >= best))
                continue;
            auto r = direct_trace::search(*branches[item.second], x, best, stats);
            if (r.cost < INF && r.cost <= best && (!found || r.cost < best)) {
                best = r.cost;
                word = move(r.word);
                found = true;
            }
        }
        // move transfers the string buffer into the thread-local result instead of
        // copying its characters. c_str() exposes its buffer but does not transfer ownership.
        if (!found)
            throw runtime_error("No feasible conditional ancestor");
        resolved = move(word);
        *score = best;
        return resolved.c_str();
    } catch (exception &e) {
        error = e.what();
        return nullptr;
    }
}
