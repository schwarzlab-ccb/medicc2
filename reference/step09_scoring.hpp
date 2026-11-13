#pragma once
// REFERENCE for migration Step 9 (not part of REF; untracked). It is REF's scoring.hpp reduced
// to the simplest correct engine: one whole-string search per WGD choice, no speed-ups.
// Left out until later steps:
//   splitting at X + shortening repeated columns  (Step 23)
//   score cache: score_key / cached_score / cache_score, the x > y swap  (Step 24)
//   bounds and branch order: BlockBounds, string_lower, witness, branch_order, room  (Steps 24/26/29)
//   pair filter + A*: make_pair_graph, greedy_pair, heuristic, budget  (Step 27)
//
// Scalar scores. score() is join()/resolve() on two single profiles, but it keeps only the
// number (search with retain = false), not the ancestors.
//   directed = true : T(x, y), x is the parent, one edge x -> y with j WGDs (k is not used).
//   directed = false: K(x, y) = min_z T(z, x) + T(z, y), two edges z -> x (k WGDs), z -> y (j WGDs).
// WGDs are whole-genome events: one (k, j) for the whole string, each WGD costs 1, charged once.

#include <algorithm>
#include <string>

#include "types.hpp"          // Model, Message, infinity
#include "fst_operations.hpp" // linear, compatible
#include "search.hpp"         // search
#include "messages.hpp"       // make_graph

// One WGD choice: the cheapest path for exactly k (left) and j (right) WGDs, without the WGD charge.
// infinity if there is no path that costs at most `cap`.
static double fixed(Model &m, const std::string &x, const std::string &y, int k, int j, bool directed,
                    double cap) {
    Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
    auto graph = make_graph(m, a, b, k, j, directed);
    return search(*graph, cap, false).cost;
}

// The minimum over all WGD choices, plus the WGD charge k + j.
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    double best = infinity;
    for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k) // directed: no edge above x
        for (int j = 0; j <= m.max_wgds; ++j)
            // Cap = best - k - j, as in join(): a branch only matters if it can reach the best so far.
            // Without it every branch searches up to its own optimum, e.g. (2, 2) explores its huge
            // product graph even after (0, 0) found cost 0 (seconds per pair instead of milliseconds).
            best = std::min(best, k + j + fixed(m, x, y, k, j, directed, best - k - j));
    return best;
}
