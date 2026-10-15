#pragma once
// Adapted from experiment 24's exact singleton-boundary factorization.
// Included inside direct_nfa after search(). No additional graph library or
// approximate pruning is introduced. Normal joins never enter this fallback.
static DFA segment_slice(const DFA &message, const vector<vector<int>> &layers,
                         int begin, int end) {
    DFA result;
    result.len = end - begin;
    vector<int> remap(message.v.size(), -1);
    for (int level = begin; level <= end; ++level)
        for (int state : layers[level]) {
            remap[state] = result.v.size();
            result.v.emplace_back();
        }
    result.start = remap[layers[begin][0]];
    for (int level = begin; level <= end; ++level)
        for (int state : layers[level]) {
            auto &node = result.v[remap[state]];
            if (level == end) {
                node.fin = end == message.len ? message.v[state].fin : 0;
                continue;
            }
            for (int symbol = 0; symbol < 10; ++symbol)
                if (message.v[state].to[symbol] >= 0) {
                    node.to[symbol] = remap[message.v[state].to[symbol]];
                    node.w[symbol] = message.v[state].w[symbol];
                }
        }
    return result;
}

static DFA segment_concatenate(const vector<DFA> &segments, int offset) {
    DFA result;
    result.offset = offset;
    for (const auto &segment : segments) {
        int shift = result.v.size();
        if (result.start < 0) result.start = shift + segment.start;
        // All optimal segment paths have zero residual weights. Splice its
        // start into previous final states, preserving original chromosome order.
        for (int state = 0; state < shift; ++state)
            if (result.v[state].fin == 0) {
                result.v[state] = segment.v[segment.start];
                for (int &target : result.v[state].to)
                    if (target >= 0) target += shift;
            }
        for (auto node : segment.v) {
            for (int &target : node.to)
                if (target >= 0) target += shift;
            result.v.push_back(node);
        }
        result.len += segment.len;
        result.offset += segment.offset;
    }
    auto layers = levels(result);
    return minimize(move(result), layers);
}

static DFA adaptive_search(AM *model, Side &a, Side &b, const vector<int> &mask,
                           int upper, int limit, Stats &stats, const vector<int> &future,
                           std::atomic<int> *incumbent, const sibling_pattern::Bound *pattern) {
    try {
        return search(a, b, mask, upper, limit, stats, future, incumbent, pattern, true);
    } catch (const SegmentRestart &restart) {
        ++stats.segment_restarts;
        const auto &cuts = restart.cuts;
        size_t count = cuts.size() - 1;
        stats.segments = count;
        auto left_layers = levels(*a.child), right_layers = levels(*b.child);
        vector<DFA> left, right, answers(count);
        left.reserve(count); right.reserve(count);
        for (size_t i = 0; i < count; ++i) {
            left.push_back(segment_slice(*a.child, left_layers, cuts[i], cuts[i + 1]));
            right.push_back(segment_slice(*b.child, right_layers, cuts[i], cuts[i + 1]));
        }
        vector<unique_ptr<Side>> xs, ys;
        vector<vector<int>> masks, futures;
        vector<unique_ptr<sibling_pattern::Bound>> patterns;
        vector<int> lower, order;
        int remaining_lower = 0;
        for (size_t i = 0; i < count; ++i) {
            xs.push_back(make_unique<Side>(left[i], a.parts, a.wg, a.wg ? 2 : 1));
            ys.push_back(make_unique<Side>(right[i], b.parts, b.wg, b.wg ? 2 : 1));
            masks.emplace_back(mask.begin() + cuts[i], mask.begin() + cuts[i + 1]);
            auto lm = lazy_cascade::layout(left[i]), rm = lazy_cascade::layout(right[i]);
            futures.push_back(chromosome_bound::future(lm, rm, 3 * a.wg + b.wg));
            patterns.push_back(a.wg <= 1 && b.wg <= 1
                ? make_unique<sibling_pattern::Bound>(model->ad, model->a, lm, rm, a.wg, b.wg)
                : nullptr);
            auto &x = *xs.back(); auto &y = *ys.back();
            int bound = max(x.h(x.start()) + y.h(y.start()),
                x.chrom_h(x.start()) + y.chrom_h(y.start()) + futures.back()[0]);
            if (patterns.back()) bound = max(bound, patterns.back()->at(0, x, x.start(), y, y.start()));
            lower.push_back(bound);
            remaining_lower = min(INF, remaining_lower + bound);
            order.push_back(i);
        }
        // Cheap factors tighten the budget before expensive ones; reconstruction
        // order is restored by storing each answer at its original index.
        stable_sort(order.begin(), order.end(), [&](int i, int j) {
            return left[i].len < left[j].len;
        });
        int wgd_cost = a.wg + b.wg;
        int offset = a.child->offset + b.child->offset + wgd_cost;
        int solved = 0;
        DFA empty; empty.len = a.child->len;
        for (int i : order) {
            upper = min(upper, incumbent->load(std::memory_order_relaxed));
            if (remaining_lower >= INF || offset + solved + remaining_lower > upper) return empty;
            remaining_lower -= lower[i];
            int room = upper - offset - solved - remaining_lower;
            auto &x = *xs[i]; auto &y = *ys[i];
            room = feasible_path::cost(x, y, masks[i], room);
            // Keep the original WGD counts for bound indexing. search() adds
            // their cost locally, so subtract it again from each segment result;
            // offset above charges it ONCE for the complete genome.
            auto result = search(x, y, masks[i], room + wgd_cost, limit, stats,
                                 futures[i], nullptr, patterns[i].get(), false);
            if (result.start < 0) return empty;
            result.offset -= wgd_cost;
            solved += result.offset;
            answers[i] = move(result);
        }
        int total = offset + solved, previous = incumbent->load(std::memory_order_relaxed);
        if (total > previous) return empty;
        while (total < previous && !incumbent->compare_exchange_weak(previous, total)) {}
        auto started = SearchClock::now();
        auto result = segment_concatenate(answers, offset);
        stats.output_seconds += elapsed(started);
        return result;
    }
}
