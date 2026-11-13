#pragma once

#include <algorithm>
#include <string>

#include "types.hpp"
#include "fst_operations.hpp"
#include "search.hpp"
#include "messages.hpp"

// Scalar scores: score() is join()/resolve() on two single profiles, keeping only the cost.
//   directed = true : T(x, y), x is the parent (one edge, j WGDs; k is not used).
//   directed = false: K(x, y) = min_z T(z, x) + T(z, y) (two edges, k and j WGDs).
// One (k, j) for the whole genome; each WGD costs 1, charged once.

// One WGD choice: the cheapest path for exactly k (left) and j (right) WGDs, without the WGD charge.
static double fixed(Model &m, const std::string &x, const std::string &y, int k, int j, bool directed) {
    Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
    auto graph = make_graph(m, a, b, k, j, directed);
    return search(*graph, infinity, false).cost;
}

// The minimum over all WGD choices, plus the WGD charge k + j
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    double best = infinity;
    for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k)
        for (int j = 0; j <= m.max_wgds; ++j)
            best = std::min(best, k + j + fixed(m, x, y, k, j, directed));
    return best;
}