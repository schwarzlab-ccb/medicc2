#pragma once
// REVIEW START: a small, relaxed version of the sibling search supplies a
// lower bound. Grouping primitive states adds possible histories; it never
// removes a real one. Each projected arc keeps the cheapest original weight.
// Child-DFA correlations and nonnegative child weights are also relaxed.
namespace sibling_pattern {
constexpr int classes = 4, side_states = 16, states = 256, unreachable = 100000000;

struct Bound {
    // Each byte stores a class number 0..3, not a pointer or a native state ID.
    vector<unsigned char> first_group, last_group;
    vector<int> values;
    int length, left_count, right_count, states;

    // A flat vector stores arcs[source][input][output][destination]. The final
    // dimension has 16 states for an abstract AD,A side. All are integer indices.
    static int arc_index(int s, int c, int o, int t) {
        return ((s * 10 + c) * 10 + o) * 16 + t;
    }

    Bound(const Part &ad, const Part &a, const vector<int> &lm, const vector<int> &rm,
          int left_wgd, int right_wgd)
        : length(lm.size()), left_count(left_wgd ? 16 : 4),
          right_count(right_wgd ? 16 : 4), states(left_count * right_count) {
        if (left_wgd < 0 || left_wgd > 1 || right_wgd < 0 || right_wgd > 1)
            throw runtime_error("sibling pattern supports zero/one WGD branches only");

        // Bit 0 means a free decrease is possible; bit 1 means a free increase
        // is possible. For AD, compare output with capped doubling of the input.
        // For A, compare with identity. This classification affects strength,
        // not safety: any partition is valid when all original arcs are projected.
        auto classify = [](const Part &p, bool doubled) {
            vector<unsigned char> groups(p.fin.size(), 0);
            for (int q = 0; q < int(p.fin.size()); ++q)
                for (int c = 0; c < 9; ++c)
                    for (const auto &e : p.arcs[q][c])
                        if (e.w == 0) {
                            int identity = doubled ? min(2 * c, 8) : c;
                            groups[q] |= (e.out < identity ? 1 : 0) |
                                         (e.out > identity ? 2 : 0);
                        }
            return groups;
        };
        first_group = classify(ad, true);
        last_group = classify(a, false);

        // The same four-dimensional layout for a single four-state primitive.
        auto ai = [](int s, int c, int o, int t) {
            return ((s * 10 + c) * 10 + o) * 4 + t;
        };
        auto quotient = [&](const Part &p, const vector<unsigned char> &groups) {
            vector<int> abstract(4 * 10 * 10 * 4, INF);
            for (int q = 0; q < int(p.fin.size()); ++q)
                for (int c = 0; c < 10; ++c)
                    for (const auto &e : p.arcs[q][c]) {
                        int &v = abstract[ai(groups[q], c, e.out, groups[e.to])];
                        v = min(v, e.w);
                    }
            return abstract;
        };
        auto first_arcs = quotient(ad, first_group), last_arcs = quotient(a, last_group);

        // Join AD and A only inside this tiny relaxation. q*4+r packs their
        // two class numbers into 0..15; mid is their shared intermediate symbol.
        vector<int> arcs(16 * 10 * 10 * 16, INF);
        for (int q = 0; q < 4; ++q)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 10; ++c)
                    for (int mid = 0; mid < 10; ++mid)
                        for (int t = 0; t < 4; ++t) {
                            int first = first_arcs[ai(q, c, mid, t)];
                            if (first >= INF)
                                continue;
                            for (int o = 0; o < 10; ++o)
                                for (int u = 0; u < 4; ++u) {
                                    int last = last_arcs[ai(r, mid, o, u)];
                                    if (last >= INF)
                                        continue;
                                    int &v = arcs[arc_index(q * 4 + r, c, o, t * 4 + u)];
                                    v = min(v, first + last);
                                }
                        }

        // values[level][joint_state] is the relaxed remaining cost. Relax final
        // weights to zero; all true final weights are nonnegative.
        values.assign(size_t(length + 1) * states, INF);
        fill(values.end() - states, values.end(), 0);
        unordered_map<int, vector<int>> side_tables, joint_tables;

        // A mask records which child symbols occur anywhere at this column.
        // Side tables retain the ancestor input c, so the siblings can agree on
        // that input later. A reference (&) borrows the cached vector, avoiding a copy.
        auto side = [&](int mask, int count) -> const vector<int> & {
            int key = mask * 16 + count - 1;
            auto it = side_tables.find(key);
            if (it != side_tables.end())
                return it->second;
            vector<int> table(count * count * 10, INF);
            for (int s = 0; s < count; ++s)
                for (int t = 0; t < count; ++t)
                    for (int c = 0; c < 10; ++c)
                        for (int o = 0; o < 10; ++o)
                            if (mask & (1 << o))
                                table[(s * count + t) * 10 + c] =
                                    min(table[(s * count + t) * 10 + c],
                                        count == 16 ? arcs[arc_index(s, c, o, t)]
                                                    : last_arcs[ai(s, c, o, t)]);
            return side_tables.emplace(key, move(table)).first->second;
        };

        // Backward dynamic programming on a layered graph. Each transition
        // consumes exactly one column, so there are no cycles to handle.
        for (int lev = length - 1; lev >= 0; --lev) {
            int code = lm[lev] * 1024 + rm[lev];
            auto found = joint_tables.find(code);
            if (found == joint_tables.end()) {
                const auto &l = side(lm[lev], left_count);
                const auto &r = side(rm[lev], right_count);
                vector<int> joint(states * states, INF);
                for (int p = 0; p < left_count; ++p)
                    for (int q = 0; q < right_count; ++q)
                        for (int t = 0; t < left_count; ++t)
                            for (int u = 0; u < right_count; ++u) {
                                int v = INF;
                                // Both sides consume the SAME ancestor symbol c.
                                for (int c = 0; c < 10; ++c)
                                    v = min(v, l[(p * left_count + t) * 10 + c] +
                                                   r[(q * right_count + u) * 10 + c]);
                                joint[(p * right_count + q) * states + t * right_count + u] = v;
                            }
                // Eviction changes only runtime: a missing matrix is recomputed.
                if (joint_tables.size() >= 64)
                    joint_tables.clear();
                found = joint_tables.emplace(code, move(joint)).first;
            }
            const auto &matrix = found->second;
            // These pointers are borrowed views of two rows in values. No new
            // memory is allocated here, and values does not resize in this loop.
            const int *next = values.data() + size_t(lev + 1) * states;
            int *current = values.data() + size_t(lev) * states;
            for (int s = 0; s < states; ++s)
                for (int t = 0; t < states; ++t)
                    current[s] = min(current[s], matrix[s * states + t] + next[t]);
        }
    }

    // The native Side key stores the primitive tuple in its high 32 bits.
    // radix decodes that tuple; the child-DFA state in the low bits is ignored
    // by this relaxation. This is numeric packing, not pointer arithmetic.
    int side_state(const direct_nfa::Side &s, uint64_t key) const {
        uint32_t code = key >> 32;
        int first = (code / s.radix.front()) % s.parts.front()->fin.size();
        int last = (code / s.radix.back()) % s.parts.back()->fin.size();
        return s.wg ? first_group[first] * 4 + last_group[last] : last_group[last];
    }

    int at(int level, const direct_nfa::Side &a, uint64_t x,
           const direct_nfa::Side &b, uint64_t y) const {
        return values[size_t(level) * states + side_state(a, x) * right_count + side_state(b, y)];
    }
};
}
