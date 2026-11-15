#pragma once
// REFERENCE for migration Step 26 (not part of REF's code). Goes to medicc/fstengine/src/scoring.hpp.
// Markers: ">>> Step 26" new or changed code, ending at "<<< Step 26" (one-line changes have no end); "--- Step 26" code to remove.
// Step 24's file; score() now uses the bounds of bounds.hpp. It is REF's score() with one change:
// the cap of a chromosome is room = min(best - other, upper), not min(best - other - 1, upper - 1),
// and there is no "certify" branch. Left out until later steps:
//   pair filter + A*: make_pair_graph  (Step 27);  greedy_pair, room = ... - 1, certify  (Step 29)

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "types.hpp"
#include "fst_operations.hpp"
#include "search.hpp"
#include "messages.hpp"
#include "bounds.hpp"

// Scalar scores: score() is join()/resolve() on two single profiles, keeping only the cost.
//   directed = true : T(x, y), x is the parent (one edge, j WGDs; k is not used).
//   directed = false: K(x, y) = min_z T(z, x) + T(z, y) (two edges, k and j WGDs).
// One (k, j) for the whole genome; each WGD costs 1, charged once.

// One WGD choice: the cheapest path for exactly k (left) and j (right) WGDs, without the WGD charge.
// infinity if there is no path that costs at most `cap`.
static double fixed(Model &m, const std::string &x, const std::string &y, int k, int j, bool directed, double cap) {
    // K is symmetric: always search in the order the cache key uses.
    if (!directed && (x > y || (x == y && k > j)))
        return fixed(m, y, x, j, k, false, cap);
    const std::string key = score_key(x, y, k, j, directed);
    const auto known = m.chromosome_scores.find(key);
    if (known != m.chromosome_scores.end())
        return known->second;
    Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
    auto graph = make_graph(m, a, b, k, j, directed);
    double value = search(*graph, cap, false).cost;
    // A finite result is the exact optimum. infinity only says "more than cap": not cached.
    if (std::isfinite(value))
        cache_score(m, key, value);
    return value;
}

// The minimum over all WGD choices of (WGD charge + sum of the chromosome scores).
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    // One entry per chromosome (the X itself is dropped: A and AD restart there). Inside a
    // chromosome, a column equal to the previous one in BOTH strings is dropped too: the
    // events that change the previous column can always cover it at no extra cost.
    std::vector<std::pair<std::string, std::string>> chromosomes(1);
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] == 'X') {
            chromosomes.emplace_back();
            continue;
        }
        auto &s = chromosomes.back();
        if (s.first.empty() || s.first.back() != x[i] || s.second.back() != y[i]) {
            s.first += x[i];
            s.second += y[i];
        }
    }
    // Shortest first: cheap chromosomes use up the budget early and stop a hopeless (k, j) sooner.
    std::stable_sort(chromosomes.begin(), chromosomes.end(),
                     [](const auto &a, const auto &b) { return a.first.size() < b.first.size(); });
    // Per chromosome and WGD choice (branch = 3k + j): a lower bound and the cost of one real
    // history (upper_bound_with_wgd). Both equal: the exact score, cached without a search.
    // >>> Step 26: change 4 lines (was 1)
    std::vector<ChromosomeBounds> bounds(chromosomes.size());
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        bounds[i].upper.fill(infinity);
        auto &[a, b] = chromosomes[i];
    // <<< Step 26
        for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k)
            for (int j = 0; j <= m.max_wgds; ++j) {
                // >>> Step 26: change 28 lines (was 4)
                int branch = 3 * k + j;
                double known = cached_score(m, a, b, k, j, directed);
                if (known >= 0) {
                    bounds[i].lower[branch] = bounds[i].upper[branch] = known;
                    continue;
                }
                double lower = directed ? 0 : lower_bound_with_wgd(a, b, branch),
                       upper = upper_bound_with_wgd(a, b, k, j, directed);
                bounds[i].lower[branch] = lower;
                bounds[i].upper[branch] = upper;
                if (lower == upper)
                    cache_score(m, score_key(a, b, k, j, directed), lower);
            }
    }
    // best starts at the cheapest upper bound; branches from the lowest lower bound.
    auto [branches, best] = branch_order(m, bounds, directed);
    for (auto [lower, branch] : branches) {
        if (lower >= best) // a scalar needs a strictly better branch
            continue;
        // bound: a lower bound of this branch that becomes exact chromosome by chromosome
        double bound = lower;
        int k = branch / 3, j = branch % 3;
        for (size_t i = 0; i < chromosomes.size(); ++i) {
            if (bound >= best)
                break;
            // other = the rest of this branch: chromosomes done + lower bounds of those to come.
            double before = bounds[i].lower[branch], other = bound - before;
            // This chromosome may use at most best - other. Its upper bound is upper, so a
            // search with cap upper always finds the exact score.
            double room = std::min(best - other, bounds[i].upper[branch]);
            double value = fixed(m, chromosomes[i].first, chromosomes[i].second, k, j, directed, room);
            if (!std::isfinite(value)) {
                bound = infinity;
                // <<< Step 26
                break;
            }
            // >>> Step 26: change 3 lines (was 1)
            bound += value - before;
        }
        best = std::min(best, bound);
            // <<< Step 26
    }
    return best;
}
