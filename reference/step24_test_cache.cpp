// REFERENCE for migration Step 24 (not part of REF's code). Goes to tests/migration/test_24_cache.cpp.
// A search stopped by a small cap must not be cached; a finished one must.
// Run: ./test_fstengine.sh 24
#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    void *handle = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    assert(handle);
    auto &m = *static_cast<Model *>(handle);
    // K("12", "21") with no WGD costs 2 (z = "11", one gain on each side), so cap 0 finds nothing ...
    assert(!std::isfinite(fixed(m, "12", "21", 0, 0, false, 0)));
    assert(m.chromosome_scores.empty()); // ... and nothing is remembered
    assert(fixed(m, "12", "21", 0, 0, false, infinity) == 2);
    assert(m.chromosome_scores.size() == 1);
    // The same pair the other way round is the same key: no second entry.
    assert(fixed(m, "21", "12", 0, 0, false, infinity) == 2);
    assert(m.chromosome_scores.size() == 1);
    // Now the cache answers even a small cap with the exact value: the caller compares it.
    assert(fixed(m, "12", "21", 0, 0, false, 0) == 2);
    assert(score(m, "12", "21", false) == 2);
    fstengine_clear(handle);
    assert(m.chromosome_scores.empty());
    assert(score(m, "12", "21", false) == 2);
    am_free(handle);
}
