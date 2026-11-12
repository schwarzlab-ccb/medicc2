// Step 7 check: search() against OpenFST's own algorithms on random graphs.
//   cost  <-> fst::ShortestDistance (no cap, no laziness)
//   tight <-> fst::Prune(f, 0): every path within 0 of the shortest
// Graphs have shared states, epsilon arcs and cycles (incl. zero-cost cycles).
#include <cassert>
#include <iostream>
#include <random>
#include "medicc/fstengine/src/search.hpp"

// Same set of paths? Encode label pairs, drop weights, then determinize + minimize both.
static bool same_paths(const Vec &a, const Vec &b) {
    fst::EncodeMapper<Arc> mapper(fst::kEncodeLabels);
    auto clean = [&](Vec f) {
        fst::ArcMap(&f, fst::RmWeightMapper<Arc>());   // every weight -> One (cost 0)
        fst::Encode(&f, &mapper);                     // transducer -> acceptor of label pairs
        fst::RmEpsilon(&f);
        Vec d;
        fst::Determinize(f, &d);
        fst::Minimize(&d);
        return d;
    };
    return fst::Equivalent(clean(a), clean(b));
}

static Vec random_graph(std::mt19937 &rng) {
    Vec g;
    int n = 2 + rng() % 7;
    for (int i = 0; i < n; ++i) g.AddState();
    g.SetStart(0);
    for (int i = 0; i < n; ++i) {
        if (rng() % 3 == 0) g.SetFinal(i, W(rng() % 3));        // some finals, cost 0-2
        int arcs = rng() % 4;
        for (int k = 0; k < arcs; ++k) {
            int label = (rng() % 5 == 0) ? 0 : 1 + rng() % 3;   // sometimes epsilon
            g.AddArc(i, Arc(label, label, W(rng() % 3), rng() % n));   // any target: cycles allowed
        }
    }
    return g;
}

int main() {
    std::mt19937 rng(7);
    int checked = 0, unreachable = 0;
    for (int round = 0; round < 3000; ++round) {
        Vec g = random_graph(rng);

        // reference cost: OpenFST's shortest distance from the start to any final
        std::vector<W> beta;
        fst::ShortestDistance(g, &beta, true);
        double best = beta.empty() ? infinity : beta[0].Value();
        if (!std::isfinite(best)) {                               // no final reachable
            assert(!std::isfinite(search(g, infinity, true).cost));
            ++unreachable;
            continue;
        }
        // reference tight set: OpenFST's Prune with threshold 0
        Vec pruned;
        fst::Prune(g, &pruned, W::One());

        for (double cap : {best - 1, best, best + 1, infinity}) {
            auto r = search(g, cap, true);
            if (cap < best) {
                assert(!std::isfinite(r.cost));
                continue;
            }
            assert(r.cost == best);                                // cost, also when cap == best
            assert(same_paths(r.tight, pruned));                   // exactly the optimal paths
            assert(search(g, cap, false).cost == best);            // cost-only mode agrees
        }
        auto [word, cost] = choose(g);                             // one optimal path
        assert(cost == best);
        Vec lw = linear(word);
        assert(std::isfinite(fst::ShortestDistance(*compose(lw, pruned)).Value()));
        ++checked;
    }
    std::cout << checked << " graphs checked against ShortestDistance + Prune(0), "
              << unreachable << " with no reachable final: PASS\n";
}
