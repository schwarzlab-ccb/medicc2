#include <cassert>
#include "medicc/fstengine/src/types.hpp"

int main() {
    assert(fst::Times(W(2), W(3)) == W(5));   // costs add along a path
    assert(fst::Plus(W(2), W(3)) == W(2));    // keep the cheaper path
    assert(W::Zero().Value() == infinity);    // "no path"
    assert(W::One().Value() == 0);            // "free"

    Message m;
    assert(m.fst.NumStates() == 0 && m.offset == 0 && m.length == 0);

    // Model model;
    // assert(model.max_wgds == 2);
    // assert(model.a.NumStates() == 0 && model.ad.NumStates() == 0);  // filled by load_model (Step 10)
    // assert(INF > 1000000);                    // NOT_FINAL in ancestral_parts.txt
}