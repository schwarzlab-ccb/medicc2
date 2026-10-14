#pragma once
#include "automata_types.hpp"
#include <boost/graph/astar_search.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <tuple>
#include <vector>

// Adapt MEDICC's bounded implicit graph to Boost.Graph's A*. The adapter owns
// no model or search nodes. Only the current expansion's arcs are materialized.
namespace library_search {
struct Edge {
    int from, to, weight, label = -1;
    bool operator==(const Edge &other) const {
        return from == other.from && to == other.to && weight == other.weight;
    }
    bool operator!=(const Edge &other) const { return !(*this == other); }
};

// Boost uses one value type for distances, edge weights, heuristics and ranks.
// Distances/edges have f == g and depth == 0. A heuristic adds the remaining
// cost, optionally prefers deeper states, and supplies the stable state-ID tie.
// The rightmost nonnegative ID wins under addition; -1 is its identity.
struct Rank {
    int f = INF, g = INF, id = -1, depth = 0;
    Rank() = default;
    Rank(int value) : f(value), g(value) {}
    Rank(int f, int g, int id, int depth = 0)
        : f(f), g(g), id(id), depth(depth) {}
};
struct Compare {
    bool operator()(Rank a, Rank b) const {
        return std::tie(a.f, a.depth, a.g, a.id) < std::tie(b.f, b.depth, b.g, b.id);
    }
};
struct Combine {
    Rank operator()(Rank a, Rank b) const {
        return {a.f + b.f, a.g + b.g, b.id < 0 ? a.id : b.id, a.depth + b.depth};
    }
};

// Borrow the vector object, not its buffer: adding implicit states may move
// its elements. Each property-map access obtains the current element again.
template<class Nodes> struct DistanceMap {
    using key_type = int;
    using value_type = Rank;
    using reference = Rank;
    using category = boost::read_write_property_map_tag;
    Nodes *nodes;
    friend Rank get(DistanceMap map, int state) {
        int g = (*map.nodes)[state].g;
        return {g, g, state};
    }
    friend void put(DistanceMap map, int state, Rank distance) {
        (*map.nodes)[state].g = distance.f;
    }
};

// Visitors use this private exception for the model's proven stopping bound.
// Other exceptions propagate to the existing C ABI error handling.
struct Stop {};
struct IgnoreRelax {
    void operator()(Edge) const {}
};
template<class Expand, class Relax> struct Graph {
    using vertex_descriptor = int;
    using edge_descriptor = Edge;
    using out_edge_iterator = std::vector<Edge>::const_iterator;
    using directed_category = boost::directed_tag;
    using edge_parallel_category = boost::allow_parallel_edge_tag;
    using traversal_category = boost::incidence_graph_tag;
    using degree_size_type = size_t;
    using vertices_size_type = size_t;
    using edges_size_type = size_t;
    using vertex_iterator = void;
    using edge_iterator = void;
    static int null_vertex() { return -1; }
    Expand expand;
    Relax relax;
    mutable std::vector<Edge> edges;
};
template<class E, class R> auto out_edges(int, const Graph<E, R> &graph) {
    return std::make_pair(graph.edges.cbegin(), graph.edges.cend());
}
template<class E, class R> size_t out_degree(int, const Graph<E, R> &graph) {
    return graph.edges.size();
}
template<class E, class R> int source(Edge edge, const Graph<E, R> &) {
    return edge.from;
}
template<class E, class R> int target(Edge edge, const Graph<E, R> &) {
    return edge.to;
}
struct Visitor : boost::default_astar_visitor {
    template<class G> void edge_relaxed(Edge edge, const G &graph) {
        graph.relax(edge);
    }
    template<class G> void examine_vertex(int state, const G &graph) {
        graph.edges.clear();
        graph.expand(state, graph.edges);
    }
};

template<class Nodes, class Expand, class H, class Relax = IgnoreRelax>
static void run(int root, Nodes &nodes, Expand expand, H heuristic, Relax relax = {}) {
    DistanceMap<Nodes> distance{&nodes};
    Graph<Expand, Relax> graph{expand, relax, {}};
    boost::vector_property_map<Rank> rank;
    boost::vector_property_map<boost::default_color_type> color;
    put(rank, root, Combine{}(get(distance, root), heuristic(root)));
    auto weight = boost::make_function_property_map<Edge>([](Edge edge) {
        return Rank{edge.weight, edge.weight, edge.to};
    });
    try {
        boost::astar_search_no_init(
            graph, root, heuristic, Visitor{}, boost::dummy_property_map{},
            rank, distance, weight, color, boost::typed_identity_property_map<int>{},
            Compare{}, Combine{}, Rank{INF}, Rank{0});
    } catch (const Stop &) {
        // The caller has already recorded its optimum or exhausted its bound.
    }
}
} // namespace library_search
