// REFERENCE for migration Step 29 (not part of REF's code). Goes to tests/migration/test_29_certify.cpp.
// The plan's test: a capped search that finds nothing is not cached; a certified upper bound is.
// Run: ./test_fstengine.sh 29
#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    void *handle = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    auto &m = *static_cast<Model *>(handle);
    // K("22", "11") with no WGD = 1 (one loss of both columns), so cap 0 finds nothing ...
    assert(!std::isfinite(fixed(m, "22", "11", 0, 0, false, 0)));
    assert(m.chromosome_scores.empty()); // ... and that is not cached
    assert(score(m, "22", "11", false) == 1);
    // score() drops the repeated column, so the chromosome is "2" / "1". Its upper bound (one loss)
    // costs 1 and nothing cheaper exists: the answer is cached (certified, or lower == upper).
    assert(m.chromosome_scores.count(score_key("2", "1", 0, 0, false)) == 1);
    assert(m.chromosome_scores.at(score_key("2", "1", 0, 0, false)) == 1);
    // Certify: for "22" / "31" with one WGD on the left, the bounds do not meet (lower 1, upper 2).
    // A search for anything cheaper than the upper bound finds nothing, so the upper bound is exact.
    assert(lower_bound_with_wgd("22", "31", 3 * 1 + 0) == 1 && upper_bound_with_wgd("22", "31", 1, 0, false) == 2);
    assert(!std::isfinite(fixed(m, "22", "31", 1, 0, false, 2 - 1)));
    assert(fixed(m, "22", "31", 1, 0, false, infinity) == 2);
    am_free(handle);
}
