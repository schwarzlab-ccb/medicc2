// REFERENCE for migration Step 26 (not part of REF's code). Goes to tests/migration/test_26_bounds.cpp.
// The plan's test, plus: the bounds really bound the exact score.
// Run: ./test_fstengine.sh 26
#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    assert(doubled("5X3", 1) == "8X6"); // capped: 10 -> 8
    assert(changes("22", "11") == 1);   // one loss over both columns
    assert(changes("2X2", "1X1") == 2); // an event never crosses X
    assert(!std::isfinite(changes("0", "1")));

    void *handle = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    auto &m = *static_cast<Model *>(handle);
    const char *pairs[][2] = {{"12", "21"}, {"1", "3"}, {"2", "4"}, {"0123", "4321"}, {"22", "11"}};
    for (auto &p : pairs)
        for (int k = 0; k <= 2; ++k)
            for (int j = 0; j <= 2; ++j) {
                double exact = fixed(m, p[0], p[1], k, j, false, infinity);
                assert(string_lower(p[0], p[1], 3 * k + j) <= exact); // lower bound
                assert(exact <= witness(p[0], p[1], k, j, false));    // upper bound
            }
    am_free(handle);
}
