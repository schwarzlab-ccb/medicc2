#pragma once
// REFERENCE for migration Step 27 (not part of REF; untracked). Goes to medicc/fstengine/src/astar.hpp (new).
// REF's astar.hpp without Matcher / Filter / StateTable / tracked_compose (they are for the event
// histories, Phase 6). branch_graph / make_graph are your functions from messages.hpp, moved here.

// Scalar pairs use cached suffix potentials for A*. Ancestor joins and directed
// searches use the same bounded OpenFST search without a separate heuristic FST:
// constructing that relaxation cost more than it saved on the measured cases.
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include <fst/project.h>

#include "types.hpp"
#include "fst_operations.hpp" // compose

// h(state): a guess of the cost still to go from a state, never too high. The guesses are
// cached per state id (NaN = not computed yet); OpenFST asks for the same state many times.
struct Heuristic {
    std::function<double(int)> evaluate;
    mutable std::vector<double> cache;
    double operator()(int state) const {
        if (state >= int(cache.size()))
            cache.resize(state + 1, std::numeric_limits<double>::quiet_NaN());
        if (std::isnan(cache[state]))
            cache[state] = evaluate(state);
        return cache[state];
    }
};

// Every parent that can reach a word of 'child' with exactly k WGDs:
// A o child, then AD o ... k times. Lazy. With `project`, keep only the parent side.
static Ptr branch_graph(Model &m, const Message &child, int k, bool project) {
    auto graph = compose(m.a, child.fst);
    for (int i = 0; i < k; ++i)
        graph = compose(m.ad, *graph);
    if (project)
        graph = std::make_unique<fst::ProjectFst<Arc>>(*graph, fst::ProjectType::INPUT);
    return graph;
}

// directed = false (join): ancestors of left with k WGDs  o  ancestors of right with j WGDs.
//   Its paths are the common ancestors z, cost T(z, left) + T(z, right) without the WGD charge.
// directed = true (resolve): parent word  o  (parent -> child edge with j WGDs, ending in the
//   child message). Input side = parent, output side = child. k is not used.
static Ptr make_graph(Model &m, const Message &left, const Message &right, int k, int j, bool directed) {
    auto r = branch_graph(m, right, j, !directed);
    auto l = directed ? Ptr(new Vec(left.fst)) : branch_graph(m, left, k, true);
    return compose(*l, *r); // ComposeFst keeps its own copies of l and r
}
