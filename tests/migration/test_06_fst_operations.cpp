#include <cassert>
#include "medicc/fstengine/src/fst_operations.hpp"

int main() {
    assert(label('0') == 1 && label('X') == 10 && symbol(9) == '8');
    Vec f = linear("1X2");
    assert(f.NumStates() == 4);
    assert(fst::ShortestDistance(f).Value() == 0);
    bool rejected = false;
    try { compatible("1X2", "12X"); } catch (const std::runtime_error &) { rejected = true; }
    assert(rejected);
}