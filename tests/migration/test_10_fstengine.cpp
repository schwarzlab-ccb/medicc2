#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    void *model = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    assert(model);
    assert(fstengine_score(model, "1X1", "2X2", 1) == 1);          // one WGD for the whole genome
    assert(!std::isfinite(fstengine_score(model, "0", "1", 1)));
    assert(fstengine_score(model, "0", "1", 0) == 1);
    void *a = am_string("1"), *b = am_string("2");
    void *joined = fstengine_join(model, a, b);
    assert(joined && am_offset(joined) == 1);
    am_drop(joined); am_drop(a); am_drop(b); am_free(model);
}