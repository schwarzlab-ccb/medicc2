#pragma once
#include <array>
#include <vector>
using std::array; using std::vector;
constexpr int INF = 100000000;
struct Node {
    array<int, 10> to, w;
    int fin = INF;
    Node() {
        to.fill(-1);
        w.fill(0);
    }
};
struct DFA {
    vector<Node> v;
    int start = -1, offset = 0, len = 0;
};

namespace direct_nfa {
struct Tight {array<vector<int>,10> to;bool final=false;};
}
