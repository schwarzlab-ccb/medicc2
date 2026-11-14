#pragma once
// REFERENCE for migration Step 23 (not part of REF; untracked). Goes to medicc/fstengine/src/scoring.hpp.
// Your Step 9 file; score() now splits the strings at X and adds chromosome scores.
// Left out until later steps:
//   score cache: score_key / cached_score / cache_score, the x > y swap in fixed()  (Step 24)
//   bounds and branch order: BlockBounds, string_lower, witness, branch_order  (Step 26)
//   pair filter + A*: make_pair_graph  (Step 27);  greedy_pair, room = ... - 1  (Step 29)

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "types.hpp"
#include "fst_operations.hpp"
#include "search.hpp"
#include "messages.hpp"

// Scalar scores: score() is join()/resolve() on two single profiles, keeping only the cost.
//   directed = true : T(x, y), x is the parent (one edge, j WGDs; k is not used).
//   directed = false: K(x, y) = min_z T(z, x) + T(z, y) (two edges, k and j WGDs).
// One (k, j) for the whole genome; each WGD costs 1, charged once.

// One WGD choice: the cheapest path for exactly k (left) and j (right) WGDs, without the WGD charge.
// infinity if there is no path that costs at most `cap`.
static double fixed(Model &m, const std::string &x, const std::string &y, int k, int j, bool directed, double cap) {
    Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
    auto graph = make_graph(m, a, b, k, j, directed);
    return search(*graph, cap, false).cost;
}

// The minimum over all WGD choices of (WGD charge + sum of the chromosome scores).
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    // One segment per chromosome (the X itself is dropped: A and AD restart there). Inside a
    // chromosome, a column equal to the previous one in BOTH strings is dropped too: the
    // events that change the previous column can always cover it at no extra cost.
    std::vector<std::pair<std::string, std::string>> segments(1);
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] == 'X') {
            segments.emplace_back();
            continue;
        }
        auto &s = segments.back();
        if (s.first.empty() || s.first.back() != x[i] || s.second.back() != y[i]) {
            s.first += x[i];
            s.second += y[i];
        }
    }
    // Shortest first: cheap segments use up the budget early and stop a hopeless (k, j) sooner.
    std::stable_sort(segments.begin(), segments.end(),
                     [](const auto &a, const auto &b) { return a.first.size() < b.first.size(); });
    double best = infinity;
    for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k) // directed: no edge above x
        for (int j = 0; j <= m.max_wgds; ++j) {
            double total = k + j;
            for (auto &[a, b] : segments) {
                // Every segment costs >= 0: this one may use at most best - total.
                total += fixed(m, a, b, k, j, directed, best - total);
                if (!std::isfinite(total))
                    break;
            }
            best = std::min(best, total);
        }
    return best;
}
