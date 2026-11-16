#pragma once
// A message is the set of optimal profiles of one tree node: an acceptor whose arc and final
// weights are all zero, plus 'offset', the subtree cost that every one of those profiles shares.
// join() bottom-up: the best common ancestors of two child messages (all ties kept)
// resolve() top-down: one best child profile for a fixed parent profile

// Optimization idea 1 -- Split at chromosomes
// A and AD restart at every X, so for a FIXED WGD choice (k, j) the cost of a profile
// is the sum of its chromosome costs. Two messages can be cut after an X where BOTH have
// exactly one state: then the pieces before and after are independent. Each piece is searched
// alone; the languages are concatenated back. The WGD choice is never made per piece: one (k, j)
// for the whole genome, charged once.

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
#include "bounds.hpp"

// The columns begin..end of a message as a message of its own (layers from fst_operations.hpp)
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

// Cut l and r into pairs of slices, one per chromosome. Neighbouring chromosomes whose paths
// are linked (more than one state after the X) stay in one pair.
static std::vector<std::pair<Message, Message>> split(const Message &l, const Message &r) {
    if (l.length != r.length)
        throw std::runtime_error("Length mismatch");
    auto l_layers = layers(l), r_layers = layers(r);
    std::vector<std::pair<Message, Message>> parts;
    int begin = 0;
    for (int end = 1; end <= l.length; ++end) {
        // Cut only after an X where both sides have exactly one state
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

// Bottom-up: the language of all minimum-cost common ancestors of l and r
// Tries every WGD count pair (k, j); each WGD costs 1, charged once per genome
static Message join(Model &m, const Message &l, const Message &r) {
    auto parts = split(l, r);
    auto bounds = message_bounds(m, parts, false);
    // Shortest chromosomes first: a short chromosome is searched faster. If the short ones
    // already use up the budget, that (k, j) is dropped without searching the long ones.
    std::vector<int> order(parts.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return parts[a].first.length < parts[b].first.length; });

    Vec combined; // union of the languages of every (k, j) that ties the best cost
    auto [branches, best] = branch_order(m, bounds, false);
    for (auto [lower, branch]: branches) {
        if (lower > best)
            continue;
        int k = branch / 3, j = branch % 3;
        double total = k + j, remaining = lower - total;
        std::vector<Vec> chromosomes(parts.size());
        bool feasible = true;
        for (int i: order) {
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            double room = best - total - remaining;
            // Every chromosome costs >= 0, so this one may use at most best - total
            // search() keeps states equal to the cap, so ties survive
            auto graph = make_graph(m, parts[i].first, parts[i].second, k, j, false);
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
        Vec language = std::move(chromosomes[0]); // chromosomes back in profile order
        for (size_t i = 1; i < chromosomes.size(); ++i)
            fst::Concat(&language, chromosomes[i]);
        if (total < best) {
            best = total;
            combined = std::move(language);
        } else if (total == best) {
            fst::Union(&combined, language);
        }
    }
    if (!std::isfinite(best) || combined.Start() < 0 || best + l.offset + r.offset >= INF)
        throw std::runtime_error("No finite common ancestor or cost exceeds range");
    return {canonical(std::move(combined)), int(best) + l.offset + r.offset, l.length};
}

// Top-down: one best child word for the fixed parent word, and its cost
// (edge cost + WGDs on this edge + the child's subtree cost in child.offset).
// One answer is enough here, so only a strictly better branch replaces it; on a tie the
// smaller j (fewer WGDs) stays.
static std::pair<std::string, double> resolve(Model &m, const std::string &parent, const Message &child) {
    Message p{linear(parent), 0, int(parent.size())};
    auto parts = split(p, child);
    auto bounds = message_bounds(m, parts, true);
    std::string bestword;
    bool found = false;
    auto [branches, best] = branch_order(m, bounds, true);

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
            auto result = search(*graph, best - total - remaining, true);
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