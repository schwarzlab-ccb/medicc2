#pragma once
#include <cstdint>
#include <vector>

namespace feasible_path {
// Pick one complete path using the EXISTING primitive transitions and suffix
// heuristic. For each ancestor symbol the two child edges are independent, so
// their best edges can be chosen without enumerating their Cartesian product.
// A failed greedy walk leaves the previous feasible upper bound unchanged.
// Excludes child offsets and WGD charges. Includes child arc/final weights.
// Symbol 9 is X, the chromosome separator. This is only a feasible upper
// bound; the caller still performs exact search, retaining ties when needed.
template<class Side>
int cost(Side &a, Side &b,
                          const vector<int> &mask, int upper) {
    uint64_t left=a.start(), right=b.start();
    int cost=0;
    for (int allowed : mask) {
        auto ap=a.prepare(left), bp=b.prepare(right);
        int best=INF, step=0;
        uint64_t next_left=0, next_right=0;
        for (int c=0; c<10; ++c) if (allowed & (1<<c)) {
            int room=upper-cost;
            if (ap.lower[0][c]+bp.lower[0][c]>room) continue;
            auto le=a.edges(ap,c,room-bp.lower[0][c]);
            if (le.empty()) continue;
            auto re=b.edges(bp,c,room-le.front().w-le.front().h);
            if (re.empty()) continue;
            const auto &x=le.front(); const auto &y=re.front();
            int value=x.w+y.w+x.h+y.h;
            if (value<best) {
                best=value;step=x.w+y.w;next_left=x.to;next_right=y.to;
            }
        }
        if (best>=INF) return upper;
        cost+=step;left=next_left;right=next_right;
    }
    return min(upper,cost+a.final(left)+b.final(right));
}
}
