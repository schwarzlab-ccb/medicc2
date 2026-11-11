#pragma once

// Scalar pairs use cached suffix potentials for A*. Ancestor joins and directed
// searches use the same bounded OpenFST search without a separate heuristic FST:
// constructing that relaxation cost more than it saved on the measured cases.
#include <functional>

struct Heuristic {
    std::function<double(int)> evaluate;
    mutable std::vector<double> cache;
    double operator()(int state) const {
        if (state >= int(cache.size()))
            cache.resize(state + 1, std::numeric_limits<double>::quiet_NaN());
        if (std::isnan(cache[state]))
            cache[state] = evaluate(state);
        return cache[state];
    }
};

using Matcher = fst::Matcher<Fst>;
using Filter = fst::SequenceComposeFilter<Matcher>;
using StateTable = fst::GenericComposeStateTable<Arc, Filter::FilterState>;

static Ptr tracked_compose(const Fst &a, const Fst &b, StateTable **table) {
    auto right = sorted(b);
    fst::ComposeFstOptions<Arc> opts(fst::CacheOptions(true, 4 << 20));
    opts.state_table = new StateTable(a, *right);
    *table = opts.state_table; // Owned by ComposeFst, including its shared copies.
    return std::make_unique<fst::ComposeFst<Arc>>(a, *right, opts);
}

static Ptr branch_graph(Model &m, const Message &child, int k, bool project) {
    auto graph = compose(m.a, child.fst);
    for (int i = 0; i < k; ++i)
        graph = compose(m.ad, *graph);
    if (project)
        graph = std::make_unique<fst::ProjectFst<Arc>>(*graph, fst::ProjectType::INPUT);
    return graph;
}

static Ptr make_graph(Model &m, const Message &left, const Message &right, int k, int j,
                      bool directed) {
    auto r = branch_graph(m, right, j, !directed);
    auto l = directed ? Ptr(new Vec(left.fst)) : branch_graph(m, left, k, true);
    return compose(*l, *r);
}
