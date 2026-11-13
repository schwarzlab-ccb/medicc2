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
};
// struct SuffixData;
// struct Model {
//     Vec a, ad, a_even, loss, gain, doubling;
//     int max_wgds = 2;
//     std::unordered_map<std::string, std::shared_ptr<const SuffixData>> suffix_cache;
//     std::unordered_map<std::string, double> chromosome_scores;
//     size_t sore_key_bytes = 0;
// };