// REFERENCE for migration Step 27 (not part of REF; untracked). Goes to tests/migration/test_27_astar.cpp.
// The plan's test, plus: the pair graph (A* + filter) gives the same answers as plain Dijkstra.
// Run: ./test_fstengine.sh 27
#include <cassert>
#include "medicc/fstengine/src/fstengine.cpp"

int main() {
    Heuristic h;
    int calls = 0;
    h.evaluate = [&](int s) { ++calls; return double(s); };
    assert(h(2) == 2 && h(2) == 2 && calls == 1); // guesses are cached

    void *handle = fstengine_load("medicc/objects/ancestral_parts.txt", 0);
    auto &m = *static_cast<Model *>(handle);
    const char *pairs[][2] = {{"12", "21"}, {"1", "3"}, {"2", "4"}, {"0123", "4321"}, {"1", "8"}, {"0", "1"}};
    for (auto &p : pairs)
        for (int k = 0; k <= 2; ++k)
            for (int j = 0; j <= 2; ++j) {
                Message a{linear(p[0]), 0, int(std::string(p[0]).size())},
                    b{linear(p[1]), 0, int(std::string(p[1]).size())};
                double dijkstra = search(*make_graph(m, a, b, k, j, false), infinity, false).cost;
                auto graph = make_pair_graph(m, p[0], p[1], k, j, infinity);
                double astar =
                    search(*graph.exact, infinity, false, graph.heuristic.get(), graph.budget.get()).cost;
                assert(astar == dijkstra);
            }
    am_free(handle);
}
