// REFERENCE for migration Step 27 (not part of REF's code). Goes to medicc/fstengine/src/types.hpp.
// Markers: ">>> Step 27" new or changed code, ending at "<<< Step 27" (one-line changes have no end); "--- Step 27" code to remove.
// Step 24's file; Model gets a_even and suffix_cache (and SuffixData is declared, it lives in
// pair_filter.hpp). Left out until later steps: score_key_bytes (Step 31); loss, gain,
// doubling (Phase 6, histories).
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
// >>> Step 27: add 1 line
struct SuffixData; // pair_filter.hpp
struct Model {
    Vec a, ad;
    // a_even: A with only even input copy numbers (and X). After a WGD every copy number is
    // even (2c, capped at 8), so the A after the last WGD can be a_even (built in fstengine.cpp).
    // >>> Step 27: add 1 line
    Vec a_even;
    int max_wgds = 2;
    // A o (one CN string), with its forward / backward distances, per string ("A"/"E" + string).
    // Built by suffix_data() (pair_filter.hpp); large, so the scheduler clears it after each row.
    // >>> Step 27: add 1 line
    std::unordered_map<std::string, std::shared_ptr<const SuffixData>> suffix_cache;
    // Exact chromosome scores, key = score_key(x, y, k, j, directed) (bounds.hpp). Only finished
    // searches are stored; a search stopped by its cap proves nothing. Cleared by fstengine_clear.
    std::unordered_map<std::string, double> chromosome_scores;
};
