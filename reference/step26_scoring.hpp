#pragma once
// REFERENCE for migration Step 26 (not part of REF; untracked). Goes to medicc/fstengine/src/scoring.hpp.
// Step 24's file; score() now uses the bounds of bounds.hpp. It is REF's score() with one change:
// the cap of a segment is room = min(best - other, upper), not min(best - other - 1, upper - 1),
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
    // Per segment and WGD choice (branch = 3k + j): a lower bound and the cost of one real
    // history (witness). Both equal: the exact score, cached without a search.
    std::vector<BlockBounds> bounds(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        bounds[i].upper.fill(infinity);
        auto &[a, b] = segments[i];
        for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k)
            for (int j = 0; j <= m.max_wgds; ++j) {
                int branch = 3 * k + j;
                double known = cached_score(m, a, b, k, j, directed);
                if (known >= 0) {
                    bounds[i].lower[branch] = bounds[i].upper[branch] = known;
                    continue;
                }
                double lower = directed ? 0 : string_lower(a, b, branch),
                       upper = witness(a, b, k, j, directed);
                bounds[i].lower[branch] = lower;
                bounds[i].upper[branch] = upper;
                if (lower == upper)
                    cache_score(m, score_key(a, b, k, j, directed), lower);
            }
    }
    // best starts at the cheapest witness; branches from the lowest lower bound.
    double best = infinity;
    auto branches = branch_order(m, bounds, directed, best);
    for (auto [lower, branch] : branches) {
        if (lower >= best) // a scalar needs a strictly better branch
            continue;
        // bound: a lower bound of this branch that becomes exact segment by segment
        double bound = lower;
        int k = branch / 3, j = branch % 3;
        for (size_t i = 0; i < segments.size(); ++i) {
            if (bound >= best)
                break;
            // other = the rest of this branch: segments done + lower bounds of those to come.
            double before = bounds[i].lower[branch], other = bound - before;
            // This segment may use at most best - other. Its witness costs upper, so a
            // search with cap upper always finds the exact score.
            double room = std::min(best - other, bounds[i].upper[branch]);
            double value = fixed(m, segments[i].first, segments[i].second, k, j, directed, room);
            if (!std::isfinite(value)) {
                bound = infinity;
                break;
            }
            bound += value - before;
        }
        best = std::min(best, bound);
    }
    return best;
}
