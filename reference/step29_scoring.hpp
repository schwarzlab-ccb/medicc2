#pragma once
// REFERENCE for migration Step 29 (not part of REF's code). Goes to medicc/fstengine/src/scoring.hpp.
// Markers: ">>> Step 29" new or changed code, ending at "<<< Step 29" (one-line changes have no end); "--- Step 29" code to remove.
// Step 27's file + greedy_pair() and "search only for something strictly cheaper". From here on
// the code is final: SAME as REF's scoring.hpp apart from comments and #includes.

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "types.hpp"
#include "fst_operations.hpp"
#include "astar.hpp"       // make_graph
#include "pair_filter.hpp" // make_pair_graph
#include "search.hpp"
#include "messages.hpp"
#include "bounds.hpp"

// Scalar scores: score() is join()/resolve() on two single profiles, keeping only the cost.
//   directed = true : T(x, y), x is the parent (one edge, j WGDs; k is not used).
//   directed = false: K(x, y) = min_z T(z, x) + T(z, y) (two edges, k and j WGDs).
// One (k, j) for the whole genome; each WGD costs 1, charged once.

// A complete path supplies a feasible upper bound, never an approximate answer.
// >>> Step 29: add 24 lines
static double greedy_pair(const PairGraph &graph) {
    // Keep active_g unset: greedy prefixes are not necessarily shortest paths.
    double cost = 0;
    for (int s = graph.exact->Start(); s >= 0;) {
        double final = graph.exact->Final(s).Value();
        if (std::isfinite(final))
            return cost + final <= graph.budget->cap ? cost + final : infinity;
        double best = infinity, weight = 0;
        int next = -1;
        for (fst::ArcIterator<Fst> it(*graph.exact, s); !it.Done(); it.Next()) {
            const auto &a = it.Value();
            double value = a.weight.Value() + (*graph.heuristic)(a.nextstate);
            if (value < best) {
                best = value;
                weight = a.weight.Value();
                next = a.nextstate;
            }
        }
        if (next < 0 || cost + best > graph.budget->cap)
            return infinity;
        cost += weight;
        s = next;
    }
    return infinity;
}
// <<< Step 29
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
    double value;
    if (!directed) {
        // Pairs: A* on the filtered pair graph (pair_filter.hpp). The filter reads the cap from
        // graph.budget, which search() lowers whenever it finds a cheaper answer.
        auto graph = make_pair_graph(m, x, y, k, j, cap);
        // Costs are whole numbers: with a feasible path of cost f, only a path of cost <= f - 1
        // can beat it, so search with that smaller cap. The answer is the cheaper of the two.
        // >>> Step 29: change 5 lines (was 1)
        double feasible = greedy_pair(graph);
        cap = std::min(cap, feasible - 1);
        graph.budget->cap = cap;
        value = std::min(feasible, search(*graph.exact, cap, false, graph.heuristic.get(),
                                          graph.budget.get()).cost);
        // <<< Step 29
    } else {
        // Directed scores stay plain bounded Dijkstra, like join() and resolve().
        Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
        auto graph = make_graph(m, a, b, k, j, directed);
        value = search(*graph, cap, false).cost;
    }
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
    std::vector<ChromosomeBounds> bounds(chromosomes.size());
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        bounds[i].upper.fill(infinity);
        auto &[a, b] = chromosomes[i];
        for (int k = 0; k <= (directed ? 0 : m.max_wgds); ++k)
            for (int j = 0; j <= m.max_wgds; ++j) {
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
            // Integral costs: only a strict improvement is needed for a scalar. This chromosome
            // only matters if it costs <= best - other - 1, and its upper bound is already upper,
            // so ask only for <= upper - 1.
            // >>> Step 29: change 1 line (was 1)
            double room = std::min(best - other - 1, bounds[i].upper[branch] - 1);
            double value = fixed(m, chromosomes[i].first, chromosomes[i].second, k, j, directed, room);
            // Nothing cheaper than the upper bound: the upper bound IS the exact score (certified).
            // >>> Step 29: add 5 lines
            if (!std::isfinite(value) && room == bounds[i].upper[branch] - 1) {
                value = bounds[i].upper[branch];
                cache_score(m, score_key(chromosomes[i].first, chromosomes[i].second, k, j, directed),
                            value);
            }
            // <<< Step 29
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
