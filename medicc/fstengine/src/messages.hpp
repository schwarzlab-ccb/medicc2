#pragma once
// A message is the set of optimal profiles of one tree node: an acceptor whose arc and final
// weights are all zero, plus 'offset', the subtree cost that every one of those profiles shares.
// join() bottom-up: the best common ancestors of two child messages (all ties kept)
// resolve() top-down: one best child profile for a fixed parent profile

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <fst/project.h>
#include <fst/union.h>

#include "types.hpp"
#include "fst_operations.hpp"
#include "search.hpp"

// Every parent that can reach a word of 'child' with exactly k WGDs:
// A o child, then AD o ... k times. Lazy. With `project`, keep only the parent side.
static Ptr branch_graph(Model &m, const Message &child, int k, bool project) {
    auto graph = compose(m.a, child.fst);
    for (int i = 0; i < k; ++i)
        graph = compose(m.ad, *graph);
    if (project)
        graph = std::make_unique<fst::ProjectFst<Arc>>(*graph, fst::ProjectType::INPUT);
    return graph;
}

// directed = false (join): ancestors of left with k WGDs  o  ancestors of right with j WGDs.
//   Its paths are the common ancestors z, cost T(z, left) + T(z, right) without the WGD charge.
// directed = true (resolve): parent word  o  (parent -> child edge with j WGDs, ending in the
//   child message). Input side = parent, output side = child. k is not used.
static Ptr make_graph(Model &m, const Message &left, const Message &right, int k, int j, bool directed) {
    auto r = branch_graph(m, right, j, !directed);
    auto l = directed ? Ptr(new Vec(left.fst)) : branch_graph(m, left, k, true);
    return compose(*l, *r); // ComposeFst keeps its own copies of l and r
}

// Bottom-up: the language of all minimum-cost common ancestors of l and r
// Tries every WGD count pair (k, j); each WGD costs 1, charged once per genome
static Message join(Model &m, const Message &l, const Message &r) {
    if (l.length != r.length)
        throw std::runtime_error("Length mismatch");
    double best = infinity;
    Vec combined; // union of the tight graphs of every (k, j) that ties the best cost
    for (int k = 0; k <= m.max_wgds; ++k) {
        for (int j = 0; j <= m.max_wgds; ++j) {
            double wgds = k + j;
            auto graph = make_graph(m, l, r, k, j, false);
            // Cap: a branch can only matter if it can tie the best so far. search() keeps
            // states that equal the cap, so ties survive. infinity also works.
            auto result = search(*graph, best - wgds, true);
            if (!std::isfinite(result.cost))
                continue;
            double total = wgds + result.cost;
            if (total < best) {
                combined.DeleteStates();
                best = total;
            }
            if (total == best) {
                fst::Union(&combined, canonical(std::move(result.tight)));
            }
        }
    }
    if (!std::isfinite(best) || combined.Start() < 0 || best + l.offset + r.offset >= INF)
        throw std::runtime_error("No finite common ancestor or cost exceeds range");
    return {canonical(std::move(combined)), best + l.offset + r.offset, l.length};
}

// Top-down: one best child word for the fixed parent word, and its cost
// (edge cost + WGDs on this edge + the child's subtree cost in child.offset).
// One answer is enough here, so only a strictly better branch replaces it.
static std::pair<std::string, double> resolve(Model &m, const std::string &parent, const Message &child) {
    if (int(parent.size()) != child.length)
        throw std::runtime_error("Length mismatch");
    Message p{linear(parent), 0, int(parent.size())};
    std::string bestword;
    double best = infinity;
    bool found = false;
    for (int j = 0; j <= m.max_wgds; ++j) {
        auto graph = make_graph(m, p, child, 0, j, true);
        auto result = search(*graph, best - j, true);
        if (!std::isfinite(result.cost))
            continue;
        double total = j + result.cost;
        if (!found || total < best) {
            best = total;
            bestword = choose(result.tight, true).first;
            found = true;
        }
    }
    if (!found)
        return {"", infinity};
    return {bestword, best + child.offset};
}