#pragma once

// Scalar scores. score() splits a pair of CN strings at X, shortens repeated
// columns per chromosome, and for every WGD branch combination sums exact
// per-chromosome scores (fixed(), cached per handle) in bound order. WGDs are
// shared across chromosomes and charged once. directed=true gives the branch
// score parent -> child; otherwise the symmetric pair distance via a common
// ancestor.

static double fixed(Model &m, const std::string &x, const std::string &y, int k, int j,
                    bool directed, double cap) {
    if (!directed && (x > y || (x == y && k > j)))
        return fixed(m, y, x, j, k, false, cap);
    const std::string key = score_key(x, y, k, j, directed);
    const auto known = m.chromosome_scores.find(key);
    if (known != m.chromosome_scores.end())
        return known->second;
    double value;
    if (!directed) {
        auto graph = make_pair_graph(m, x, y, k, j, cap);
        value = search(*graph.exact, cap, false, graph.heuristic.get(), graph.budget.get()).cost;
    } else {
        Message a{linear(x), 0, int(x.size())}, b{linear(y), 0, int(y.size())};
        auto graph = make_graph(m, a, b, k, j, directed);
        value = search(*graph.exact, cap, false, graph.heuristic.get()).cost;
    }
    // A bounded failure proves only value > cap, never global infeasibility.
    if (std::isfinite(value))
        m.chromosome_scores.emplace(key, value);
    return value;
}
static double score(Model &m, const std::string &x, const std::string &y, bool directed) {
    compatible(x, y);
    // Fixed WGD counts are shared across chromosomes and charged once.
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
    std::stable_sort(segments.begin(), segments.end(),
                     [](const auto &a, const auto &b) { return a.first.size() < b.first.size(); });
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
                    m.chromosome_scores.emplace(score_key(a, b, k, j, directed), lower);
            }
    }
    double best = infinity;
    auto branches = branch_order(m, bounds, directed, best);
    for (auto [lower, branch] : branches) {
        if (lower >= best)
            continue;
        double bound = lower;
        int k = branch / 3, j = branch % 3;
        for (size_t i = 0; i < segments.size(); ++i) {
            if (bound >= best)
                break;
            double before = bounds[i].lower[branch], other = bound - before;
            // Integral costs: only a strict improvement is needed for a scalar.
            double room = std::min(best - other - 1, bounds[i].upper[branch]);
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
