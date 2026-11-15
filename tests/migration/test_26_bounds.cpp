#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    assert(doubled("5X3", 1) == "8X6"); // capped: 10 -> 8
    assert(upper_bound_no_wgd("22", "11") == 1);   // one loss over both columns
    assert(upper_bound_no_wgd("2X2", "1X1") == 2); // an event never crosses X
    assert(!std::isfinite(upper_bound_no_wgd("0", "1")));

    void *handle = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    auto &m = *static_cast<Model *>(handle);
    const char *pairs[][2] = {{"12", "21"}, {"1", "3"}, {"2", "4"}, {"0123", "4321"}, {"22", "11"}};
    for (auto &p : pairs)
        for (int k = 0; k <= 2; ++k)
            for (int j = 0; j <= 2; ++j) {
                double exact = fixed(m, p[0], p[1], k, j, false, infinity);
                assert(lower_bound_with_wgd(p[0], p[1], 3 * k + j) <= exact); // lower bound
                assert(exact <= upper_bound_with_wgd(p[0], p[1], k, j, false));    // upper bound
            }
    am_free(handle);
}
