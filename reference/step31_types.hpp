// REFERENCE for migration Step 31 (not part of REF; untracked). Goes to medicc/fstengine/src/types.hpp.
// Step 27's file + score_key_bytes. Code now SAME as REF's types.hpp except loss, gain,
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
struct SuffixData; // pair_filter.hpp
struct Model {
    Vec a, ad;
    // a_even: A with only even input copy numbers (and X). After a WGD every copy number is
    // even (2c, capped at 8), so the A after the last WGD can be a_even (built in fstengine.cpp).
    Vec a_even;
    int max_wgds = 2;
    // A o (one CN string), with its forward / backward distances, per string ("A"/"E" + string).
    // Built by suffix_data() (pair_filter.hpp); large, so the scheduler clears it after each row.
    std::unordered_map<std::string, std::shared_ptr<const SuffixData>> suffix_cache;
    // Exact chromosome scores, key = score_key(x, y, k, j, directed) (bounds.hpp). Only finished
    // searches are stored; a search stopped by its cap proves nothing. Bounded by cache_score; kept across rows.
    std::unordered_map<std::string, double> chromosome_scores;
    size_t score_key_bytes = 0; // total length of the keys in chromosome_scores (cache_score)
};
