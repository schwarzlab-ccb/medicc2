#pragma once

// Symbol mapping, linear acceptors for CN strings, the layer decomposition of
// fixed-length messages and delayed (lazy) OpenFST composition helpers.

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
        throw std::runtime_error("Invalid label");
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
        label(x[i]);
        label(y[i]);
        if ((x[i] == 'X') != (y[i] == 'X'))
            throw std::runtime_error("Separator mismatch");
    }
}
struct Layers {
    std::vector<std::vector<int>> states;
    std::vector<int> depth;
    std::vector<bool> after_x;
};
static Layers layers(const Message &m) {
    Layers out;
    out.states.resize(m.length + 1);
    out.depth.assign(m.fst.NumStates(), -1);
    out.after_x.assign(m.length + 1, true);
    if (m.fst.Start() < 0)
        throw std::runtime_error("Empty message");
    out.states[0].push_back(m.fst.Start());
    out.depth[m.fst.Start()] = 0;
    out.after_x[0] = false;
    for (int d = 0; d <= m.length; ++d)
        for (int state : out.states[d]) {
            if (m.fst.Final(state) != W::Zero() && d != m.length)
                throw std::runtime_error("Message is not fixed length");
            for (fst::ArcIterator<Vec> it(m.fst, state); !it.Done(); it.Next()) {
                const auto &a = it.Value();
                if (d == m.length || a.ilabel != a.olabel || a.ilabel < 1 || a.ilabel > 10)
                    throw std::runtime_error("Invalid message arc");
                if (a.ilabel != 10)
                    out.after_x[d + 1] = false;
                if (out.depth[a.nextstate] < 0) {
                    out.depth[a.nextstate] = d + 1;
                    out.states[d + 1].push_back(a.nextstate);
                } else if (out.depth[a.nextstate] != d + 1)
                    throw std::runtime_error("Message has inconsistent path lengths");
            }
        }
    return out;
}

static Ptr sorted(const Fst &f) {
    return std::make_unique<fst::ArcSortFst<Arc, fst::ILabelCompare<Arc>>>(
        f, fst::ILabelCompare<Arc>());
}
static Ptr compose(const Fst &a, const Fst &b) {
    // Sort only arcs of states that are requested. Never enumerate either operand.
    auto right = sorted(b);
    fst::ComposeFstOptions<Arc> opts(fst::CacheOptions(true, 4 << 20));
    return std::make_unique<fst::ComposeFst<Arc>>(a, *right, opts);
}
