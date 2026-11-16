#pragma once
// REFERENCE for migration Step 27, OPTIONAL: Heuristic written as a `class` instead of a `struct`.
// For reading only; step27_astar.hpp (the struct) stays the version to type. Both behave the same.
//
// What changes:
//   struct: `evaluate` and `cache` are public; the builder writes  h.evaluate = [...](int s) {...};
//   class : both are private; the formula goes in through the constructor, and only h(s) and
//           size() are open to other code. Nobody can replace the formula or edit the cache.
// What has to change elsewhere if you use this version:
//   pair_filter.hpp, make_pair_graph():
//       auto heuristic = std::make_shared<Heuristic>();
//       heuristic->evaluate = [table, l = left.suffix, r = right.suffix](int state) {...};
//     becomes
//       auto heuristic = std::make_shared<Heuristic>([table, l = left.suffix, r = right.suffix](int state) {...});
//   test_27_astar.cpp:
//       Heuristic h;  h.evaluate = [&](int s) {...};
//     becomes
//       Heuristic h([&](int s) { ++calls; return double(s); });
//   search.hpp needs no change: it only calls (*heuristic)(s).
#include <cmath>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

// h(state): a guess of the cost still to go from a state, never too high. The guesses are
// cached per state id (NaN = not computed yet); OpenFST asks for the same state many times.
class Heuristic {
  public:
    // explicit: a lambda does not silently turn into a Heuristic, you write Heuristic(...) yourself
    explicit Heuristic(std::function<double(int)> evaluate) : evaluate_(std::move(evaluate)) {}

    // h(state): computed on the first call, read from the cache after that
    double operator()(int state) const {
        if (state >= int(cache_.size()))
            cache_.resize(state + 1, std::numeric_limits<double>::quiet_NaN());
        if (std::isnan(cache_[state]))
            cache_[state] = evaluate_(state);
        return cache_[state];
    }

    // How many state ids have a slot in the cache (for tests; not needed by the engine)
    int size() const { return int(cache_.size()); }

  private:
    // Trailing _ is a common C++ habit for private members (REF uses it in CostComposeFilter)
    std::function<double(int)> evaluate_; // the formula, fixed at construction
    mutable std::vector<double> cache_;   // mutable: operator() is const but fills the cache
};
