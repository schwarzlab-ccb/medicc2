#include <cassert>
#include <vector>
#include "medicc/fstengine/src/fst_operations.hpp"

// layers() on a small hand-drawn message (labels are arbitrary; 10 = X):
//
//   s0 -1-> s1 -2-> s2 -X-> s3 -3-> s4 -4-> s5 -X-> s6 -5-> s7 (final)
//            \             ^  \             ^         \       ^
//             3-> s8 --X--/    9-> s9 --10-/           --10--/
//
//   layer:    0  1  2      3  4      5  6  7
//   after_x:  F  F  F      T  F      F  T  F

static bool throws(const Message &m) {
    try { layers(m); } catch (const std::runtime_error &) { return true; }
    return false;
}

int main() {
    const int X = 10;
    Vec f;
    for (int s = 0; s < 10; ++s)
        f.AddState();
    f.SetStart(0);
    f.SetFinal(7, W::One());
    auto arc = [&](int from, int l, int to) { f.AddArc(from, Arc(l, l, W::One(), to)); };
    arc(0, 1, 1); arc(1, 2, 2); arc(2, X, 3); arc(3, 3, 4); arc(4, 4, 5); arc(5, X, 6); arc(6, 5, 7);
    arc(1, 3, 8); arc(8, X, 3);
    arc(3, 9, 9); arc(9, 10, 5);
    arc(6, 10, 7);
    Message m{f, 0, 7};

    Layers ls = layers(m);
    // states: layer -> states (order = order found)
    std::vector<std::vector<int>> states = {{0}, {1}, {2, 8}, {3}, {4, 9}, {5}, {6}, {7}};
    assert(ls.states == states);
    // depth: state -> layer, the reverse of states
    std::vector<int> depth = {0, 1, 2, 3, 4, 5, 6, 7, 2, 4};
    assert(ls.depth == depth);
    // after_x: one entry per LAYER (length + 1 = 8), not per state
    std::vector<bool> after_x = {false, false, false, true, false, false, true, false};
    assert(ls.after_x == after_x);
    for (int s = 0; s < f.NumStates(); ++s) // the per-state view is after_x[depth[s]]
        assert(ls.after_x[ls.depth[s]] == (s == 3 || s == 6));

    // A final state before the last layer: not fixed length.
    Message early = m;
    early.fst.SetFinal(4, W::One());
    assert(throws(early));
    // W::Zero() (= infinity) means "not final": setting it changes nothing.
    Message zero = m;
    zero.fst.SetFinal(4, W::Zero());
    assert(!throws(zero));
    // Two paths of different length into the same state.
    Message uneven = m;
    uneven.fst.AddArc(0, Arc(1, 1, W::One(), 2));
    assert(throws(uneven));
}
