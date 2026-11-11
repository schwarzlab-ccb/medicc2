#pragma once

// A* potentials from OpenFST state abstractions. Each primitive is mapped to a
// small state quotient (Pattern) and each child message to its layer graph; the
// relaxed product of these quotients is materialised and its backward shortest
// distances give a consistent heuristic for the delayed exact product, because
// every concrete arc has an abstract counterpart of no greater weight.
#include <functional>

struct Heuristic {
    std::function<double(int)> evaluate;
    std::function<int(int)> project;
    std::vector<W> backward;
    int states = 0;
    mutable std::vector<double> cache;
    double operator()(int state) const {
        if (state >= int(cache.size()))
            cache.resize(state + 1, std::numeric_limits<double>::quiet_NaN());
        if (!std::isnan(cache[state]))
            return cache[state];
        if (evaluate)
            return cache[state] = evaluate(state);
        int abstract = project(state);
        // Every concrete reachable state maps to an abstract reachable state.
        if (abstract < 0 || abstract >= states)
            throw std::runtime_error("A* abstraction missed a reachable state");
        double value = abstract < int(backward.size()) ? backward[abstract].Value() : infinity;
        return cache[state] = value;
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

static Ptr sum_arcs(const Fst &f) {
    // Remove redundant parallel arcs on the small relaxation, per requested state.
    return std::make_unique<fst::StateMapFst<Arc, Arc, fst::ArcSumMapper<Arc>>>(
        f, fst::ArcSumMapper<Arc>(f));
}

struct Lattice {
    Ptr exact, relaxed;
    std::function<int(int)> project;
};

static Lattice product(Lattice left, Lattice right) {
    StateTable *actual_table, *relaxed_table;
    auto exact = tracked_compose(*left.exact, *right.exact, &actual_table);
    auto relaxed = tracked_compose(*left.relaxed, *right.relaxed, &relaxed_table);
    auto map = [actual_table, relaxed_table, l = left.project, r = right.project](int s) {
        // Copy the tuple: recursive projection can grow other state tables.
        auto t = actual_table->Tuple(s);
        return relaxed_table->FindState(
            typename StateTable::StateTuple(l(t.StateId1()), r(t.StateId2()), t.GetFilterState()));
    };
    return {std::move(exact), sum_arcs(*relaxed), std::move(map)};
}

static Lattice child_lattice(const Message &message) {
    auto ls = layers(message);
    Vec relaxed;
    for (int d = 0; d <= message.length; ++d)
        relaxed.AddState();
    relaxed.SetStart(0);
    relaxed.SetFinal(message.length, W::One());
    for (int d = 0; d < message.length; ++d) {
        int mask = 0;
        for (int s : ls.states[d])
            for (fst::ArcIterator<Vec> it(message.fst, s); !it.Done(); it.Next())
                mask |= 1 << it.Value().ilabel;
        for (int l = 1; l <= 10; ++l)
            if (mask & (1 << l))
                relaxed.AddArc(d, Arc(l, l, W::One(), d + 1));
    }
    return {Ptr(new Vec(message.fst)), Ptr(new Vec(std::move(relaxed))),
            [depth = std::move(ls.depth)](int s) { return depth.at(s); }};
}

static Pattern quotient(const Vec &primitive, bool doubled, int classes) {
    Pattern pattern;
    pattern.group.resize(primitive.NumStates());
    for (int s = 0; s < primitive.NumStates(); ++s) {
        int bits = 0;
        for (fst::ArcIterator<Vec> it(primitive, s); !it.Done(); it.Next()) {
            auto a = it.Value();
            if (a.ilabel == 10 || a.weight != W::One())
                continue;
            int cn = a.ilabel - 1, identity = doubled ? std::min(2 * cn, 8) : cn;
            bits |= (a.olabel - 1 < identity ? 1 : 0) | (a.olabel - 1 > identity ? 2 : 0);
        }
        // Four classes retain independent free-loss/free-gain flags. With two
        // WGDs, coarsen to two classes to keep the pattern graph inexpensive.
        pattern.group[s] = classes == 4 ? bits : ((bits & 2) ? 1 : 0);
    }
    for (int s = 0; s < classes; ++s)
        pattern.fst.AddState();
    pattern.fst.SetStart(pattern.group[primitive.Start()]);
    for (int s = 0; s < primitive.NumStates(); ++s) {
        int source = pattern.group[s];
        pattern.fst.SetFinal(source, fst::Plus(pattern.fst.Final(source), primitive.Final(s)));
        for (fst::ArcIterator<Vec> it(primitive, s); !it.Done(); it.Next()) {
            auto a = it.Value();
            a.nextstate = pattern.group[a.nextstate];
            pattern.fst.AddArc(source, a);
        }
    }
    fst::ArcSumMapper<Arc> mapper(pattern.fst);
    Vec compact;
    fst::StateMap(pattern.fst, &compact, mapper);
    pattern.fst = std::move(compact);
    fst::ArcSort(&pattern.fst, fst::ILabelCompare<Arc>());
    return pattern;
}

static Lattice primitive_lattice(const Vec &primitive, const Pattern &pattern) {
    return {Ptr(new Vec(primitive)), Ptr(new Vec(pattern.fst)),
            [groups = pattern.group](int s) { return groups.at(s); }};
}

static Lattice branch_lattice(Model &m, const Message &child, int k, int classes, bool project) {
    const Pattern &a = classes == 4 ? m.a4 : m.a2;
    const Pattern &ad = classes == 4 ? m.ad4 : m.ad2;
    auto p = product(primitive_lattice(m.a, a), child_lattice(child));
    for (int i = 0; i < k; ++i)
        p = product(primitive_lattice(m.ad, ad), std::move(p));
    if (project) {
        p.exact = std::make_unique<fst::ProjectFst<Arc>>(*p.exact, fst::ProjectType::INPUT);
        p.relaxed = std::make_unique<fst::ProjectFst<Arc>>(*p.relaxed, fst::ProjectType::INPUT);
    }
    return p;
}

struct Graph {
    Ptr exact, relaxed; // Keep abstract state tables alive for the projection.
    std::shared_ptr<Heuristic> heuristic;
};

static Graph make_graph(Model &m, const Message &left, const Message &right, int k, int j,
                        bool directed) {
    int classes = (k == 2 || j == 2) ? 2 : 4;
    auto r = branch_lattice(m, right, j, classes, !directed);
    auto l = directed ? child_lattice(left) : branch_lattice(m, left, k, classes, true);
    auto graph = product(std::move(l), std::move(r));
    auto heuristic = std::make_shared<Heuristic>();
    heuristic->project = std::move(graph.project);
    // This materializes only the small state quotient constrained by the child
    // columns. The exact graph and full P remain delayed throughout.
    Vec small(*graph.relaxed);
    heuristic->states = small.NumStates();
    fst::ShortestDistance(small, &heuristic->backward, true);
    check(small);
    if (heuristic->backward.size() == 1 && !heuristic->backward[0].Member())
        throw std::runtime_error("Pattern shortest distance failed");
    return {std::move(graph.exact), std::move(graph.relaxed), std::move(heuristic)};
}
