#pragma once

// Lower and upper bounds that order and prune the exact searches. Upper bounds
// are concrete feasible histories (upper_bound_with_wgd); lower bounds relax the capped
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
// Bound both entry count and retained string bytes. Clearing only evicts exact
// memoized answers; it never changes a search bound or caches a cutoff failure.
static void cache_score(Model &m, std::string key, double value) {
    constexpr size_t max_entries = 65536, max_key_bytes = 8 << 20;
    if (key.size() > max_key_bytes)
        return;
    if (m.chromosome_scores.size() >= max_entries ||
        m.score_key_bytes + key.size() > max_key_bytes) {
        m.chromosome_scores.clear();
        m.score_key_bytes = 0;
    }
    auto inserted = m.chromosome_scores.emplace(std::move(key), value);
    if (inserted.second)
        m.score_key_bytes += inserted.first->first.size();
}

static double cached_score(Model &m, const std::string &x, const std::string &y, int k, int j,
                           bool directed) {
    auto found = m.chromosome_scores.find(score_key(x, y, k, j, directed));
    return found == m.chromosome_scores.end() ? -1 : found->second;
}
// Upper bound for ONE A stage when there is NO WGD (A = L_LOH o G: losses, then gains).
// It is the cost of one concrete history in which every column only loses or only gains.
// Losses and gains are counted apart: pl / pg = the loss / gain depth of the previous column.
// Pay 1 for every NEW event (a column deeper than the one before it); ending events is free.
// X ends all events. Infinity: different lengths, X in only one string, or 0 -> positive.
// Example (gains above the line, losses below; one # = one event over that column):
//
//   column    1    2    3    4    X    6    7    8
//   source    3    3    2    1    X    1    2    2
//   target    1    2    2    3    X    2    4    3
//
//   gain +2                  #    |         #
//   gain +1                  #    |    #    #    #
//          -----------------------+---------------
//   loss -1   #    #              |
//   loss -2   #                   |
//
//   new       2    0    0    2    |    1    1    0    total 6
//
// Not exact: a running loss cannot carry over a 0 -> 0 column, so 202 -> 000 gives 4
// (the model needs 2: two losses over all three columns, the middle 0 stays 0).
// upper_bound_with_wgd() adds the WGDs (doubled) and the ancestors around this function.
static double upper_bound_no_wgd(const std::string &source, const std::string &target) {
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
// The smallest copy number that is >= cn after k doublings: cn / 2^k, rounded up.
// Example: cn = 5, k = 1 -> 3 (3 doubled is 6 >= 5; 2 doubled is only 4).
static int smallest_before_doublings(int cn, int k) {
    int factor = 1; // 2^k
    for (int i = 0; i < k; ++i)
        factor *= 2;
    return (cn + factor - 1) / factor; // integer division, rounded up
}
// The upper bound WITH WGDs: a concrete feasible history with identity A stages, the
// mandatory capped doublings, then upper_bound_no_wgd in the last A. No zero-to-positive gain.
static double upper_bound_with_wgd(const std::string &x, const std::string &y, int k, int j, bool directed) {
    if (directed)
        return upper_bound_no_wgd(doubled(x, j), y);
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
        best = std::min(best, upper_bound_no_wgd(doubled(ancestor, k), x) + upper_bound_no_wgd(doubled(ancestor, j), y));
    }
    // One more guess, not flat: in each column, the smallest ancestor copy number that is still
    // >= x after k doublings and >= y after j doublings. Then neither branch needs a gain.
    for (size_t i = 0; i < x.size(); ++i)
        if (x[i] != 'X') {
            int need_x = smallest_before_doublings(x[i] - '0', k);
            int need_y = smallest_before_doublings(y[i] - '0', j);
            ancestor[i] = char('0' + std::max(need_x, need_y));
        }
    return std::min(best, upper_bound_no_wgd(doubled(ancestor, k), x) + upper_bound_no_wgd(doubled(ancestor, j), y));
}
static double lower_bound_with_wgd(const std::string &x, const std::string &y, int branch) {
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
// The bits of a column mask (column_masks): bit c = copy number c can be here.
constexpr int X_BIT = 1 << 9;      // bit 9: the column can be X
constexpr int CN_BITS = X_BIT - 1; // bits 0..8: the column can be copy number 0..8
static double lower_bound_message(const std::vector<int> &left,
                                  const std::vector<int> &right, int branch) {
    // Drop child correlations and nonnegative weights ONLY for this lower bound.
    // Native VectorFst and ShortestDistance replace the production handwritten DP.
    for (size_t i = 0; i < left.size(); ++i) {
        bool lx = left[i] & X_BIT, rx = right[i] & X_BIT;
        if (lx != rx || (lx && ((left[i] & CN_BITS) || (right[i] & CN_BITS))))
            return 0;
    }
    double total = 0;
    for (size_t first = 0; first < left.size();) {
        if (left[first] == X_BIT) {
            ++first;
            continue;
        }
        size_t end = first;
        while (end < left.size() && left[end] != X_BIT)
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
                    (directed ? 0 : lower_bound_message(lm, rm, branch)) + l.second + r.second;
                bounds.upper[branch] =
                    upper_bound_with_wgd(l.first, r.first, k, j, directed) + l.second + r.second;
            }
    }
    return out;
}
// Returns the WGD choices sorted from the lowest lower bound, and the cheapest upper
// bound: the best cost the caller starts with.
static std::pair<std::vector<std::pair<double, int>>, double>
branch_order(Model &model, const std::vector<BlockBounds> &blocks, bool directed) {
    std::vector<std::pair<double, int>> order;
    double best = infinity;
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
    return {order, best};
}
