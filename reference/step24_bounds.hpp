#pragma once
// REFERENCE for migration Step 24 (not part of REF's code). Goes to medicc/fstengine/src/bounds.hpp (new).
// The first two functions of REF's bounds.hpp: the key of the score cache and a store without a
// size limit. Left out until later steps:
//   cached_score, changes / doubled / witness, string_lower, language_lower, BlockBounds,
//   message_bounds, branch_order, #include "chromosome_constants.hpp"  (Step 26)
//   the size limit of cache_score (65536 entries, 8 MiB of keys, score_key_bytes)  (Step 31)
//
// The cache remembers exact chromosome scores: a distance matrix asks for the same
// chromosome pairs again and again (same profiles, other partner).

#include <string>
#include <utility>

#include "types.hpp"

// "K" or "D", the WGD counts, then the two strings. K(x, y) with (k, j) equals K(y, x) with
// (j, k), so undirected keys are stored in one order only: x <= y.
static std::string score_key(std::string x, std::string y, int k, int j, bool directed) {
    if (!directed && (x > y || (x == y && k > j))) {
        std::swap(x, y);
        std::swap(k, j);
    }
    return std::string(directed ? "D" : "K") + char('0' + k) + char('0' + j) + x + ":" + y;
}

// Store one EXACT score. Never store a search that stopped at its cap: "nothing within the
// cap" is not "impossible". No size limit yet; the scheduler clears the cache after each row.
static void cache_score(Model &m, std::string key, double value) {
    m.chromosome_scores.emplace(std::move(key), value);
}
