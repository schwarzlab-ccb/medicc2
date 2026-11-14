// REFERENCE for migration Step 23 (not part of REF; untracked). Goes to tests/migration/test_23_split.cpp.
// The plan's test, plus empty chromosomes ("1XX2") and a check that split blocks add up.
// Run: ./test_fstengine.sh 23
#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    Message a{linear("1X2"), 0, 3}, b{linear("2X3"), 0, 3};
    auto blocks = split(a, b);
    assert(blocks.size() == 2);
    assert(blocks[0].first.length == 2 && blocks[1].first.length == 1); // "1X" | "2"

    // Empty chromosome: a cut after each X, the middle block is just "X".
    Message c{linear("1XX2"), 0, 4}, d{linear("2XX3"), 0, 4};
    blocks = split(c, d);
    assert(blocks.size() == 3);
    assert(blocks[1].first.length == 1);

    void *model = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    assert(model);
    auto &m = *static_cast<Model *>(model);
    // An empty chromosome costs nothing, and here the best (k, j) is the same per chromosome ...
    assert(score(m, "1XX2", "2XX3", true) == score(m, "1", "2", true) + score(m, "2", "3", true));
    // ... but one WGD is shared by the whole genome: 1 -> 2 on both chromosomes costs 1, not 2.
    assert(score(m, "1X1", "2X2", true) == 1);
    // join(): the offset is the same as the scalar score of the two leaves.
    Message x{linear("12X21"), 0, 5}, y{linear("21X12"), 0, 5};
    assert(join(m, x, y).offset == score(m, "12X21", "21X12", false));
    am_free(model);
}
