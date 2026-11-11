#pragma once

// Lower and upper bounds that order and prune the exact searches. Upper bounds
// are concrete feasible histories (witness); lower bounds relax the capped
// L_LOH model at segment boundaries (tables in chromosome_constants.hpp, made
// by tools/generate_chromosome_bound.py). Branches are the 3x3 combinations of
// WGD counts on the two sides of a join or pair score.
#include "chromosome_constants.hpp"
#include <array>
#include <numeric>

static std::string score_key(std::string x, std::string y, int k, int j, bool directed) {
    if (!directed && (x > y || (x == y && k > j))) {
        std::swap(x, y);
        std::swap(k, j);
    }
    return std::string(directed ? "D" : "K") + char('0' + k) + char('0' + j) + x + ":" + y;
}
static double cached_score(Model &m, const std::string &x, const std::string &y, int k, int j,
                           bool directed) {
    auto found = m.chromosome_scores.find(score_key(x, y, k, j, directed));
    return found == m.chromosome_scores.end() ? -1 : found->second;
}
// A concrete feasible history: identity A stages, mandatory capped doublings,
// then unit interval losses/gains in the last A. No zero-to-positive gain.
static double changes(const std::string &source, const std::string &target) {
    if (source.size() != target.size())
        return infinity;
    int pl = 0, pg = 0, total = 0;
    for (size_t i = 0; i < source.size(); ++i) {
        if (source[i] == 'X' || target[i] == 'X') {
            if (source[i] != target[i])
                return infinity;
            pl = pg = 0;
            continue;
        }
        int a = source[i] - '0', b = target[i] - '0';
        if (a == 0 && b > 0)
            return infinity;
        int loss = std::max(0, a - b), gain = std::max(0, b - a);
        total += std::max(0, loss - pl) + std::max(0, gain - pg);
        pl = loss;
        pg = gain;
    }
    return total;
}
static std::string doubled(std::string word, int k) {
    for (char &c : word)
        if (c != 'X')
            c = '0' + std::min(8, (c - '0') * (1 << k));
    return word;
}
static double witness(const std::string &x, const std::string &y, int k, int j, bool directed) {
    if (directed)
        return changes(doubled(x, j), y);
    if (x.size() != y.size())
        return infinity;
    double best = infinity;
    std::string ancestor = x;
    for (int c = 1; c <= 8; ++c) {
        for (size_t i = 0; i < x.size(); ++i) {
            if ((x[i] == 'X') != (y[i] == 'X'))
                return infinity;
            ancestor[i] = x[i] == 'X' ? 'X' : char('0' + c);
        }
        best = std::min(best, changes(doubled(ancestor, k), x) + changes(doubled(ancestor, j), y));
    }
    for (size_t i = 0; i < x.size(); ++i)
        if (x[i] != 'X')
            ancestor[i] = '0' + std::max(((x[i] - '0') + (1 << k) - 1) / (1 << k),
                                         ((y[i] - '0') + (1 << j) - 1) / (1 << j));
    return std::min(best, changes(doubled(ancestor, k), x) + changes(doubled(ancestor, j), y));
}
static double string_lower(const std::string &x, const std::string &y, int branch) {
    // A single path through the boundary relaxation; no graph search needed.
    if (x.empty())
        return 0;
    int total = 0, previous = -1;
    for (size_t i = 0; i < x.size(); ++i) {
        int u = x[i] - '0', v = y[i] - '0', state = 9 * u + v;
        total += previous < 0
                     ? chromosome_bound::END[branch][state]
                     : chromosome_bound::EDGE[branch]
                                             [(previous / 9 * 9 + u) * 81 + previous % 9 * 9 + v];
        previous = state;
    }
    return (total + chromosome_bound::END[branch][previous] + 1) / 2;
}
static std::vector<int> column_masks(const Message &message) {
    auto ls = layers(message);
    std::vector<int> masks(message.length, 0);
    for (int d = 0; d < message.length; ++d)
        for (int s : ls.states[d])
            for (fst::ArcIterator<Vec> it(message.fst, s); !it.Done(); it.Next())
                masks[d] |= 1 << (it.Value().ilabel - 1);
    return masks;
}
static double language_lower(const std::vector<int> &left,
                             const std::vector<int> &right, int branch) {
    // Drop child correlations and nonnegative weights ONLY for this lower bound.
    // Native VectorFst and ShortestDistance replace the production handwritten DP.
    for (size_t i = 0; i < left.size(); ++i) {
        bool lx = left[i] & 512, rx = right[i] & 512;
        if (lx != rx || (lx && ((left[i] & 511) || (right[i] & 511))))
            return 0;
    }
    double total = 0;
    for (size_t first = 0; first < left.size();) {
        if (left[first] == 512) {
            ++first;
            continue;
        }
        size_t end = first;
        while (end < left.size() && left[end] != 512)
            ++end;
        Vec graph;
        graph.SetStart(graph.AddState());
        std::vector<std::pair<int, int>> previous;
        for (size_t i = first; i < end; ++i) {
            std::vector<std::pair<int, int>> current;
            for (int x = 0; x < 9; ++x)
                if (left[i] & (1 << x))
                    for (int y = 0; y < 9; ++y)
                        if (right[i] & (1 << y)) {
                            int state = 9 * x + y, id = graph.AddState();
                            current.emplace_back(state, id);
                            if (i == first)
                                graph.AddArc(0, Arc(state + 1, state + 1,
                                                    W(chromosome_bound::END[branch][state]), id));
                            else
                                for (auto [old, source] : previous)
                                    graph.AddArc(
                                        source,
                                        Arc(state + 1, state + 1,
                                            W(chromosome_bound::EDGE[branch]
                                                                    [(old / 9 * 9 + x) * 81 +
                                                                     old % 9 * 9 + y]),
                                            id));
                        }
            previous = std::move(current);
        }
        for (auto [state, id] : previous)
            graph.SetFinal(id, W(chromosome_bound::END[branch][state]));
        double cost = fst::ShortestDistance(graph).Value();
        check(graph);
        total += std::ceil(cost / 2);
        first = end;
    }
    return total;
}
struct BlockBounds {
    std::array<double, 9> lower{}, upper{};
};
static std::vector<BlockBounds>
message_bounds(Model &model, const std::vector<std::pair<Message, Message>> &parts, bool directed) {
    std::vector<BlockBounds> out(parts.size());
    for (size_t i = 0; i < parts.size(); ++i) {
        auto &bounds = out[i];
        bounds.upper.fill(infinity);
        auto l = choose(parts[i].first.fst), r = choose(parts[i].second.fst);
        auto lm = column_masks(parts[i].first), rm = column_masks(parts[i].second);
        for (int k = 0; k <= (directed ? 0 : model.max_wgds); ++k)
            for (int j = 0; j <= model.max_wgds; ++j) {
                int branch = 3 * k + j;
                bounds.lower[branch] =
                    (directed ? 0 : language_lower(lm, rm, branch)) + l.second + r.second;
                bounds.upper[branch] =
                    witness(l.first, r.first, k, j, directed) + l.second + r.second;
            }
    }
    return out;
}
static std::vector<std::pair<double, int>>
branch_order(Model &model, const std::vector<BlockBounds> &blocks, bool directed, double &best) {
    std::vector<std::pair<double, int>> order;
    for (int k = 0; k <= (directed ? 0 : model.max_wgds); ++k)
        for (int j = 0; j <= model.max_wgds; ++j) {
            int branch = 3 * k + j;
            double lower = k + j, upper = k + j;
            for (auto &block : blocks) {
                lower += block.lower[branch];
                upper += block.upper[branch];
            }
            if (upper < best)
                best = upper;
            order.emplace_back(lower, branch);
        }
    std::stable_sort(order.begin(), order.end(), [](const auto &a, const auto &b) {
        if (a.first != b.first)
            return a.first < b.first;
        int ac = a.second / 3 + a.second % 3, bc = b.second / 3 + b.second % 3;
        return ac == bc ? a.second < b.second : ac < bc;
    });
    return order;
}
