#pragma once

// Symbol mapping, linear acceptors for CN strings, the layer decomposition of
// fixed-length message and dealyed (lazy) OpenFST composition helpers

#include "types.hpp"

static thread_local std::string error, output;
static int label(char c) {
    if (c == 'X')
        return 10;
    if (c >= '0' && c <= '8')
        return c - '0' + 1;
    throw std::runtime_error("Invalid symbol");
}
static char symbol(int l) {
    if (l < 1 || l > 10)
        throw std::runtime_error("Invalid symbol");
    return l == 10 ? 'X' : char('0' + l - 1);
}
static void check(const Fst &f) {
    if (f.Properties(fst::kError, false))
        throw std::runtime_error("OpenFST error");
}
static Vec linear(const std::string &s) {
    Vec f;
    for (size_t i = 0; i <= s.size(); ++i)
        f.AddState();
    f.SetStart(0);
    f.SetFinal(s.size(), W::One());
    for (size_t i = 0; i < s.size(); ++i)
        f.AddArc(i, Arc(label(s[i]), label(s[i]), W::One(), i + 1));
    return f;
}
static void compatible(const std::string &x, const std::string &y) {
    if (x.size() != y.size())
        throw std::runtime_error("Length mismatch");
    for (size_t i = 0; i < x.size(); ++i) {
        label(x[i]); // check if x[i] is a valid label
        label(y[i]);
        if ((x[i] == 'X') != (y[i] == 'X'))
            throw std::runtime_error("Separator mismatch");
    }
}
//TODO: add layers struct and static layers here
static Ptr sorted(const Fst &f) {
    return std::make_unique<fst::ArcSortFst<Arc, fst::ILabelCompare<Arc>>>(
        f, fst::ILabelCompare<Arc>());
}
static Ptr compose(const Fst &a, const Fst &b) {
    //Sort only arcs of states that are requested. Never enumerate either operand
    auto right = sorted(b);
    fst::ComposeFstOptions<Arc> opts(fst::CacheOptions(true, 4 << 20));
    return std::make_unique<fst::ComposeFst<Arc>>(a, *right, opts);
}


