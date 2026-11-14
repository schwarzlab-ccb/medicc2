#pragma once
// REFERENCE for migration Step 26 (not part of REF's code). Goes to medicc/fstengine/src/messages.hpp.
// Markers: ">>> Step 26" new or changed code, ending at "<<< Step 26" (one-line changes have no end); "--- Step 26" code to remove.
// Step 23's file; join() and resolve() now use the bounds. join() / resolve() are SAME as REF's
// messages.hpp. Left out until later steps: make_graph / branch_graph move to astar.hpp (Step 27).
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
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fst/concat.h>
#include <fst/project.h>
#include <fst/union.h>

#include "types.hpp"
#include "fst_operations.hpp"
#include "search.hpp"
// >>> Step 26: add 1 line
#include "bounds.hpp"   // message_bounds, branch_order

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

// Cut two equally long messages into chromosomes (pairs of slices). Neighbouring chromosomes
// whose paths are linked (more than one state at the X between them) stay in one piece.
static std::vector<std::pair<Message, Message>> split(const Message &l, const Message &r) {
    if (l.length != r.length)
        throw std::runtime_error("Length mismatch");
    auto l_layers = layers(l), r_layers = layers(r);
    std::vector<std::pair<Message, Message>> parts;
    int begin = 0;
    for (int end = 1; end <= l.length; ++end) {
        // A unique state on each side after X proves independence; otherwise
        // keep the correlations by keeping these chromosomes together.
        if (end != l.length &&
            !(l_layers.states[end].size() == 1 && r_layers.states[end].size() == 1 &&
              l_layers.after_x[end] && r_layers.after_x[end]))
            continue;
        parts.emplace_back(slice(l, l_layers, begin, end), slice(r, r_layers, begin, end));
        begin = end;
    }
    if (l.length == 0)
        parts.emplace_back(l, r);
    return parts;
}

// Bottom-up: the language of all minimum-cost common ancestors of l and r.
// Bounds (bounds.hpp) order the WGD choices from the lowest lower bound and drop a choice, or
// stop it early, once its lower bound is STRICTLY worse than the best cost: equal choices are
// kept, because their ancestors are optimal too. `remaining` = lower bounds of the chromosomes not
// searched yet; a chromosome may use at most best - total - remaining.
static Message join(Model &m, const Message &l, const Message &r) {
    auto parts = split(l, r);
    // >>> Step 26: add 1 line
    auto bounds = message_bounds(m, parts, false);
    std::vector<int> order(parts.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return parts[a].first.length < parts[b].first.length; });
    double best = infinity;
    Vec combined;
    // >>> Step 26: change 6 lines (was 3)
    auto branches = branch_order(m, bounds, false, best);
    for (auto [lower, branch] : branches) {
        if (lower > best)
            continue; // Equality must retain all optimal words.
        int k = branch / 3, j = branch % 3;
        double total = k + j, remaining = lower - total;
    // <<< Step 26
        std::vector<Vec> chromosomes(parts.size());
        bool feasible = true;
        for (int i : order) {
            // >>> Step 26: add 6 lines
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            double room = best - total - remaining;
            // <<< Step 26
            auto graph = make_graph(m, parts[i].first, parts[i].second, k, j, false);
            // >>> Step 26: change 1 line (was 1)
            auto result = search(*graph, room, true);
            if (!std::isfinite(result.cost)) {
                feasible = false;
                break;
            }
            total += result.cost;
            chromosomes[i] = canonical(std::move(result.tight));
        }
        if (!feasible)
            continue;
        Vec language = std::move(chromosomes[0]);
        for (size_t i = 1; i < chromosomes.size(); ++i)
            fst::Concat(&language, chromosomes[i]);
        if (total < best) {
            combined.DeleteStates();
            best = total;
        }
        if (total == best)
            fst::Union(&combined, language);
    }
    if (!std::isfinite(best) || combined.Start() < 0 || best + l.offset + r.offset >= INF)
        throw std::runtime_error("No finite common ancestor or cost exceeds range");
    // >>> Step 26: change 1 line (was 1)
    return {canonical(std::move(combined)), int(best) + l.offset + r.offset, l.length};
}
// Top-down: one best child word for the fixed parent word, and its cost
// (edge cost + WGDs on this edge + the child's subtree cost in child.offset).
// One answer is enough here, so a choice that can only tie the best is skipped too.
// >>> Step 26: change 2 lines (was 1)
static std::pair<std::string, double> resolve(Model &m, const std::string &parent,
                                              const Message &child) {
// <<< Step 26
    Message p{linear(parent), 0, int(parent.size())};
    auto parts = split(p, child);
    // >>> Step 26: add 1 line
    auto bounds = message_bounds(m, parts, true);
    std::string bestword;
    double best = infinity;
    bool found = false;
    // >>> Step 26: change 5 lines (was 2)
    auto branches = branch_order(m, bounds, true, best);
    for (auto [lower, branch] : branches) {
        if (lower > best || (found && lower == best))
            continue;
        double total = branch, remaining = lower - total;
    // <<< Step 26
        std::string word;
        bool feasible = true;
        // >>> Step 26: change 9 lines (was 3)
        for (size_t i = 0; i < parts.size(); ++i) {
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            auto graph = make_graph(m, parts[i].first, parts[i].second, 0, branch, true);
            auto result =
                search(*graph, best - total - remaining, true);
        // <<< Step 26
            if (!std::isfinite(result.cost)) {
                feasible = false;
                break;
            }
            total += result.cost;
            word += choose(result.tight, true).first;
        }
        if (feasible && (!found || total < best)) {
            best = total;
            bestword = std::move(word);
            found = true;
        }
    }
    if (!found)
        return {"", infinity};
    return {bestword, best + child.offset};
}
