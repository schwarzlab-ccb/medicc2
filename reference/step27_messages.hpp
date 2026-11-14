#pragma once
// REFERENCE for migration Step 27 (not part of REF's code). Goes to medicc/fstengine/src/messages.hpp.
// Markers: ">>> Step 27" new or changed code, ending at "<<< Step 27" (one-line changes have no end); "--- Step 27" code to remove.
// Step 26's file without branch_graph / make_graph (moved to astar.hpp). From here on the code
// is final: SAME as REF's messages.hpp apart from comments and #includes.
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
// --- Step 27: remove 1 line here, starting: #include <fst/project.h>
#include <fst/union.h>

#include "types.hpp"
#include "fst_operations.hpp"
// >>> Step 27: add 1 line
#include "astar.hpp"    // make_graph
#include "search.hpp"
#include "bounds.hpp"   // message_bounds, branch_order

// The columns begin..end of a message as a message of its own (layers from fst_operations.hpp).
// --- Step 27: remove 13 lines here, starting: static Ptr branch_graph(Model &m, const Message &child, int k, bool pr
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
    if (l.length != r.length)
        throw std::runtime_error("Length mismatch");
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
// Bounds (bounds.hpp) order the WGD choices from the lowest lower bound and drop a choice, or
// stop it early, once its lower bound is STRICTLY worse than the best cost: equal choices are
// kept, because their ancestors are optimal too. `remaining` = lower bounds of the blocks not
// searched yet; a block may use at most best - total - remaining.
static Message join(Model &m, const Message &l, const Message &r) {
    auto parts = split(l, r);
    auto bounds = message_bounds(m, parts, false);
    std::vector<int> order(parts.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return parts[a].first.length < parts[b].first.length; });
    double best = infinity;
    Vec combined;
    auto branches = branch_order(m, bounds, false, best);
    for (auto [lower, branch] : branches) {
        if (lower > best)
            continue; // Equality must retain all optimal words.
        int k = branch / 3, j = branch % 3;
        double total = k + j, remaining = lower - total;
        std::vector<Vec> blocks(parts.size());
        bool feasible = true;
        for (int i : order) {
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            double room = best - total - remaining;
            auto graph = make_graph(m, parts[i].first, parts[i].second, k, j, false);
            auto result = search(*graph, room, true);
            if (!std::isfinite(result.cost)) {
                feasible = false;
                break;
            }
            total += result.cost;
            blocks[i] = canonical(std::move(result.tight));
        }
        if (!feasible)
            continue;
        Vec language = std::move(blocks[0]);
        for (size_t i = 1; i < blocks.size(); ++i)
            fst::Concat(&language, blocks[i]);
        if (total < best) {
            combined.DeleteStates();
            best = total;
        }
        if (total == best)
            fst::Union(&combined, language);
    }
    if (!std::isfinite(best) || combined.Start() < 0 || best + l.offset + r.offset >= INF)
        throw std::runtime_error("No finite common ancestor or cost exceeds range");
    return {canonical(std::move(combined)), int(best) + l.offset + r.offset, l.length};
}
// Top-down: one best child word for the fixed parent word, and its cost
// (edge cost + WGDs on this edge + the child's subtree cost in child.offset).
// One answer is enough here, so a choice that can only tie the best is skipped too.
static std::pair<std::string, double> resolve(Model &m, const std::string &parent,
                                              const Message &child) {
    Message p{linear(parent), 0, int(parent.size())};
    auto parts = split(p, child);
    auto bounds = message_bounds(m, parts, true);
    std::string bestword;
    double best = infinity;
    bool found = false;
    auto branches = branch_order(m, bounds, true, best);
    for (auto [lower, branch] : branches) {
        if (lower > best || (found && lower == best))
            continue;
        double total = branch, remaining = lower - total;
        std::string word;
        bool feasible = true;
        for (size_t i = 0; i < parts.size(); ++i) {
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            auto graph = make_graph(m, parts[i].first, parts[i].second, 0, branch, true);
            auto result =
                search(*graph, best - total - remaining, true);
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
