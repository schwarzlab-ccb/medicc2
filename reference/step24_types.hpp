// REFERENCE for migration Step 24 (not part of REF's code). Goes to medicc/fstengine/src/types.hpp.
// Markers: ">>> Step 24" new or changed code, ending at "<<< Step 24" (one-line changes have no end); "--- Step 24" code to remove.
// Your current file; Model gets the score cache `chromosome_scores` (your commented-out copy of
// REF's Model is dropped). Left out until later steps: a_even, suffix_cache (Step 27),
// score_key_bytes (Step 31); loss, gain, doubling (Phase 6, histories).
#pragma once

// Shared types of the fstengine (copy-number engine). A Message is a fixed-length
// weighted acceptor over the CN alphabet 0..8 and X;
// its integer offset carries cost that was pushed out of the acceptor. The
// Model holds the model primitives (A = L_LOH o G, AD = A o D with the capped
// doubling D) read from objects/ancestral_parts.txt, and per-handle caches.
// A Model is never shared by threads.
#include <fst/fstlib.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using W = fst::TropicalWeightTpl<double>;
using Arc = fst::ArcTpl<W>;
using Fst = fst::Fst<Arc>;
using Vec = fst::VectorFst<Arc>;
using Ptr = std::unique_ptr<Fst>;
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr int INF = 100000000;
struct Message {
    Vec fst;
    double offset = 0;
    int length = 0;
};
struct Model {
    Vec a, ad;
    int max_wgds = 2;
    // Exact chromosome scores, key = score_key(x, y, k, j, directed) (bounds.hpp). Only finished
    // searches are stored; a search stopped by its cap proves nothing. Cleared by fstengine_clear.
    // >>> Step 24: add 1 line
    std::unordered_map<std::string, double> chromosome_scores;
};
