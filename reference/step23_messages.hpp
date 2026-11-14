#pragma once
// REFERENCE for migration Step 23 (not part of REF's code). Goes to medicc/fstengine/src/messages.hpp.
// Markers: ">>> Step 23" new or changed code, ending at "<<< Step 23" (one-line changes have no end); "--- Step 23" code to remove.
// Your Step 8 file + slice() / split() (copied from REF's messages.hpp), and join() / resolve()
// search block by block. Left out until later steps:
//   bounds: message_bounds, branch_order, the `remaining` / `lower` pruning (Step 26)
//   make_graph / branch_graph move to astar.hpp (Step 27)
//
// A message is the set of optimal profiles of one tree node: an acceptor whose arc and final
// weights are all zero, plus 'offset', the subtree cost that every one of those profiles shares.
// join() bottom-up: the best common ancestors of two child messages (all ties kept)
// resolve() top-down: one best child profile for a fixed parent profile
//
// Splitting: A and AD restart at every X, so for a FIXED WGD choice (k, j) the cost of a profile
// is the sum of its chromosome costs. Two messages can be cut after an X where BOTH have exactly
// one state: then the pieces before and after are independent (no correlation crosses the cut).
// Each piece is searched alone; the languages are concatenated back. The WGD choice is never
// made per piece: one (k, j) for the whole genome, charged once.

#include <cmath>
#include <memory>
// >>> Step 23: add 1 line
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
// >>> Step 23: add 2 lines
#include <vector>

#include <fst/concat.h>
// <<< Step 23
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

// The columns begin..end of a message as a message of its own (layers from fst_operations.hpp).
// >>> Step 23: change 23 lines (was 1)
static Message slice(const Message &m, const Layers &ls, int begin, int end) {
    Message out;
    out.length = end - begin;
    std::vector<int> ids(m.fst.NumStates(), -1);
    for (int d = begin; d <= end; ++d)
        for (int s : ls.states[d])
            ids[s] = out.fst.AddState();
    out.fst.SetStart(ids[ls.states[begin].at(0)]);
    for (int d = begin; d <= end; ++d)
        for (int s : ls.states[d]) {
            if (d == end) {
                out.fst.SetFinal(ids[s], end == m.length ? m.fst.Final(s) : W::One());
                continue;
            }
            for (fst::ArcIterator<Vec> it(m.fst, s); !it.Done(); it.Next()) {
                auto a = it.Value();
                a.nextstate = ids[a.nextstate];
                out.fst.AddArc(ids[s], a);
            }
        }
    return out;
}

// Cut two equally long messages into independent blocks (pairs of slices).
static std::vector<std::pair<Message, Message>> split(const Message &l, const Message &r) {
// <<< Step 23
    if (l.length != r.length)
        throw std::runtime_error("Length mismatch");
    // >>> Step 23: add 20 lines
    auto a = layers(l), b = layers(r);
    std::vector<std::pair<Message, Message>> parts;
    int begin = 0;
    for (int end = 1; end <= l.length; ++end) {
        // A unique state on each side after X proves independence; otherwise
        // keep the correlations by leaving the block whole.
        if (end != l.length && !(a.states[end].size() == 1 && b.states[end].size() == 1 &&
                                 a.after_x[end] && b.after_x[end]))
            continue;
        parts.emplace_back(slice(l, a, begin, end), slice(r, b, begin, end));
        begin = end;
    }
    if (l.length == 0)
        parts.emplace_back(l, r);
    return parts;
}

// Bottom-up: the language of all minimum-cost common ancestors of l and r.
// Tries every WGD count pair (k, j); for each, searches the blocks one by one and adds.
static Message join(Model &m, const Message &l, const Message &r) {
    auto parts = split(l, r);
    // Shortest blocks first: a cheap block uses up the budget early and stops a hopeless
    // (k, j) sooner.
    std::vector<int> order(parts.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return parts[a].first.length < parts[b].first.length; });
    // <<< Step 23
    double best = infinity;
    Vec combined; // union of the languages of every (k, j) that ties the best cost
    // >>> Step 23: change 1 line (was 1)
    for (int k = 0; k <= m.max_wgds; ++k)
        for (int j = 0; j <= m.max_wgds; ++j) {
            // >>> Step 23: change 14 lines (was 4)
            double total = k + j;
            std::vector<Vec> blocks(parts.size());
            bool feasible = true;
            for (int i : order) {
                // Every block costs >= 0, so this block may use at most best - total.
                // search() keeps states equal to the cap, so ties survive.
                auto graph = make_graph(m, parts[i].first, parts[i].second, k, j, false);
                auto result = search(*graph, best - total, true);
                if (!std::isfinite(result.cost)) {
                    feasible = false;
                    break;
                }
                total += result.cost;
                blocks[i] = canonical(std::move(result.tight));
            }
            if (!feasible)
            // <<< Step 23
                continue;
            // >>> Step 23: change 3 lines (was 1)
            Vec language = std::move(blocks[0]); // blocks back in profile order
            for (size_t i = 1; i < blocks.size(); ++i)
                fst::Concat(&language, blocks[i]);
            // <<< Step 23
            if (total < best) {
                combined.DeleteStates();
                best = total;
            }
            // >>> Step 23: change 2 lines (was 4)
            if (total == best)
                fst::Union(&combined, language);
            // <<< Step 23
        }
    if (!std::isfinite(best) || combined.Start() < 0 || best + l.offset + r.offset >= INF)
        throw std::runtime_error("No finite common ancestor or cost exceeds range");
    return {canonical(std::move(combined)), best + l.offset + r.offset, l.length};
}

// Top-down: one best child word for the fixed parent word, and its cost
// (edge cost + WGDs on this edge + the child's subtree cost in child.offset).
// One answer is enough here, so only a strictly better branch replaces it.
static std::pair<std::string, double> resolve(Model &m, const std::string &parent, const Message &child) {
    // --- Step 23: remove 2 lines here, starting: if (int(parent.size()) != child.length)
    Message p{linear(parent), 0, int(parent.size())};
    // >>> Step 23: add 1 line
    auto parts = split(p, child);
    std::string bestword;
    double best = infinity;
    bool found = false;
    for (int j = 0; j <= m.max_wgds; ++j) {
        // >>> Step 23: change 14 lines (was 6)
        double total = j;
        std::string word;
        bool feasible = true;
        for (auto &part : parts) {
            auto graph = make_graph(m, part.first, part.second, 0, j, true);
            auto result = search(*graph, best - total, true);
            if (!std::isfinite(result.cost)) {
                feasible = false;
                break;
            }
            total += result.cost;
            word += choose(result.tight, true).first;
        }
        if (feasible && (!found || total < best)) {
        // <<< Step 23
            best = total;
            // >>> Step 23: change 1 line (was 1)
            bestword = std::move(word);
            found = true;
        }
    }
    if (!found)
        return {"", infinity};
    return {bestword, best + child.offset};
}
