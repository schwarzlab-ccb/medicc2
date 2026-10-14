// OpenFST owns the standard automata algorithms. Only this translation unit
// sees its templates; search kernels see the small value-type adapter API.
#include "openfst_language.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <fst/determinize.h>
#include <fst/push.h>
#include <fst/rmepsilon.h>
#include <fst/union.h>

using namespace std;
namespace library_language {
using Weight = fst::TropicalWeightTpl<double>;
using Arc = fst::ArcTpl<Weight>;
using Fst = fst::VectorFst<Arc>;

static int integer(Weight weight) {
    double value = weight.Value();
    if (value == numeric_limits<double>::infinity()) return INF;
    if (!isfinite(value) || value < 0 || value >= INF || floor(value) != value)
        throw runtime_error("OpenFST returned a non-integral or out-of-range cost");
    return int(value);
}

// Read the borrowed MEDICC acceptor; return a separately owned OpenFST value.
// OpenFST reserves label zero for epsilon, so MEDICC symbols 0..9 map to 1..10.
static Fst encode(const DFA &message) {
    Fst result;
    for (size_t state = 0; state < message.v.size(); ++state) result.AddState();
    if (message.start < 0) return result;
    result.SetStart(message.start);
    for (int state = 0; state < int(message.v.size()); ++state) {
        const auto &node = message.v[state];
        if (node.fin < INF) result.SetFinal(state, Weight(node.fin));
        for (int symbol = 0; symbol < 10; ++symbol) {
            if (node.to[symbol] >= 0) {
                int label = symbol + 1;
                result.AddArc(state, Arc(label, label, Weight(node.w[symbol]), node.to[symbol]));
            }
        }
    }
    return result;
}

static DFA decode(const Fst &acceptor, int length, int offset) {
    if (acceptor.Properties(fst::kError, true))
        throw runtime_error("OpenFST language operation failed");
    DFA result;
    result.len = length;
    if (acceptor.Start() == fst::kNoStateId) return result;
    result.start = acceptor.Start();
    result.offset = offset;
    result.v.resize(acceptor.NumStates());
    for (int state = 0; state < acceptor.NumStates(); ++state) {
        result.v[state].fin = integer(acceptor.Final(state));
        for (fst::ArcIterator<Fst> arcs(acceptor, state); !arcs.Done(); arcs.Next()) {
            const auto &arc = arcs.Value();
            int symbol = arc.ilabel - 1;
            if (symbol < 0 || symbol >= 10 || arc.ilabel != arc.olabel || result.v[state].to[symbol] >= 0)
                throw runtime_error("OpenFST result is not an epsilon-free deterministic acceptor");
            result.v[state].to[symbol] = arc.nextstate;
            result.v[state].w[symbol] = integer(arc.weight);
        }
    }
    return result;
}

DFA determinized(const vector<direct_nfa::Tight> &nodes, int start, int length,
                 int offset, size_t &states) {
    Fst input, result;
    for (size_t state = 0; state < nodes.size(); ++state) input.AddState();
    input.SetStart(start);
    for (int state = 0; state < int(nodes.size()); ++state) {
        if (nodes[state].final) input.SetFinal(state, Weight::One());
        for (int symbol = 0; symbol < 10; ++symbol) {
            for (int target : nodes[state].to[symbol]) {
                int label = symbol + 1;
                input.AddArc(state, Arc(label, label, Weight::One(), target));
            }
        }
    }
    // &result lets OpenFST fill this local object; no ownership is transferred.
    fst::Determinize(input, &result, fst::DeterminizeOptions<Arc>(fst::kDelta));
    states = result.NumStates();
    return decode(result, length, offset);
}

DFA united(const vector<DFA> &messages) {
    Fst combined, result;
    int length = 0;
    for (const auto &message : messages) {
        length = message.len;
        if (message.start < 0) continue;
        auto acceptor = encode(message);
        // Tropical RemoveWeight subtracts its argument: subtracting -offset
        // puts the message's stored cost back onto its accepting paths.
        fst::RemoveWeight(&acceptor, Weight(-message.offset), false);
        fst::Union(&combined, acceptor);
    }
    fst::RmEpsilon(&combined);
    fst::Determinize(combined, &result, fst::DeterminizeOptions<Arc>(fst::kDelta));
    return decode(result, length, 0);
}
}  // namespace library_language
