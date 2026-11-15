#pragma once

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
    if (std::isfinite(value))
        cache_score(m, key, value);
    return value;
}

// The minimum over all WGD choices of (WGD charge + sum of the chromosome scores).
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    // One entry per chromosome: [(chr_1_x, chr_1_y), (chr_2_x, chr_2_y), ...] (the X is dropped:
    // A and AD restart there). Inside a chromosome, a column equal to the previous one in BOTH
    // strings is dropped too: the events that change the previous column can cover it at no
    // extra cost.
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
    // Per chromosome and WGD choice (branch = 3k + j): a lower bound and the cost of one real history
    // (upper_bound_with_wgd). Both equal: the exact score, cached without a search
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

    auto [branches, best] = branch_order(m, bounds, directed);
    for (auto [lower, branch] : branches) {
        if (lower >= best) // a scaler needs a strictly better branch
            continue;
        // bound: a lower bound of this branch that becomes exact chromosome by chromosome
        double bound = lower;
        int k = branch / 3, j = branch % 3;
        for (size_t i = 0; i < chromosomes.size(); ++i) {
            if (bound >= best)
                break;
            double before = bounds[i].lower[branch], other = bound - before;
            double room = std::min(best - other, bounds[i].upper[branch]);
            double value = fixed(m, chromosomes[i].first, chromosomes[i].second, k, j, directed, room);
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
