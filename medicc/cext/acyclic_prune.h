#ifndef MEDICC_FSTLIB_EXTENSION_ACYCLIC_PRUNE_H_
#define MEDICC_FSTLIB_EXTENSION_ACYCLIC_PRUNE_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fst/connect.h>
#include <fst/expanded-fst.h>
#include <fst/fstlib.h>
#include <fst/matcher.h>
#include <fst/prune.h>
#include <fst/script/fst-class.h>

namespace medicc_fstlib_extension {

inline fst::script::MutableFstClass* prune_acyclic_std(
    const fst::script::FstClass& input_class,
    float weight_threshold, float delta) {

    using Arc = fst::StdArc;
    using StateId = Arc::StateId;
    using Weight = Arc::Weight;

    if (input_class.ArcType() != Arc::Type()) {
        throw std::invalid_argument(
            "prune_acyclic requires an OpenFST standard-arc FST");
    }

    const fst::Fst<Arc>* input = input_class.GetFst<Arc>();
    if (input == nullptr) {
    throw std::runtime_error("could not obtain the typed OpenFST input");
    }

    constexpr std::uint64_t required = fst::kAcyclic | fst::kTopSorted;
    const std::uint64_t properties = input->Properties(required, true);
    if ((properties & required) != required) {
    throw std::invalid_argument(
        "prune_acyclic requires an acyclic, topologically sorted FST");
    }

    const StateId num_states = fst::CountStates(*input);
    std::vector<Weight> distance(num_states, Weight::Zero());

    // In the tropical semiring used by StdArc:
    //
    //   fst::Plus(a, b)  means min(a, b)
    //   fst::Times(a, b) means a + b
    //   Weight::One()    is numeric cost 0
    //   Weight::Zero()   is numeric cost +infinity
    //
    // For each state s, the desired recurrence is:
    //
    //   distance[s] = min(
    //       final_weight[s],
    //       arc.weight + distance[arc.nextstate] for every outgoing arc
    //   )
    //
    // final_weight[s] represents the cost of ending a successful path at s.  An
    // ordinary non-final state has Final(s) == Weight::Zero(), so it contributes
    // no finite candidate until an outgoing arc reaches a coaccessible state.
    //
    // Because kTopSorted guarantees arc.nextstate > state, scanning state IDs in
    // descending order guarantees every distance[arc.nextstate] was computed in
    // an earlier loop iteration.  This is the piece the generic reverse shortest-
    // distance implementation cannot assume for an arbitrary FST.

    for (StateId state=num_states; state-- >0;) {
        Weight best = input->Final(state);
        for (fst::ArcIterator<fst::Fst<Arc>> arcs(*input, state); !arcs.Done(); arcs.Next()) {
            const    Arc& arc = arcs.Value();
            const Weight candidate = fst::Times(arc.weight, distance[arc.nextstate]);
            best = fst::Plus(best, candidate);
        }
        distance[state] = best;
    }

    auto output = std::make_unique<fst::VectorFst<Arc>>();


    // The constructor arguments are, in order:
    //
    // 1. pruning weight threshold;
    // 2. no maximum state-count threshold;
    // 3. accept every arc;
    // 4. our precomputed distance-to-final vector;
    // 5. numerical delta.
    //
    // Supplying argument 4 is the entire memory-saving hook.  If it were nullptr,
    // fst::Prune would invoke ShortestDistance(input, reverse=true), and that
    // generic reverse path materializes the enormous reversed VectorFst.
    fst::PruneOptions<Arc, fst::AnyArcFilter<Arc>> options(
        Weight(weight_threshold),
        fst::kNoStateId,
        fst::AnyArcFilter<Arc>(),
        &distance,
        delta
    );

    fst::Prune(*input, output.get(), options);
    return new fst::script::VectorFstClass(std::move(output));
}

namespace internal{

using FusedArc = fst::StdArc;
using FusedFst = fst::Fst<FusedArc>;
using FusedStateId = FusedArc::StateId;
using FusedWeight = FusedArc::Weight;

// A product state is just a pair of ordinary OpenFST state IDs.  Packing the
// pair into one integer makes the hash table substantially smaller than an
// unordered_map<std::pair<...>, ...> implementation.
// instead of <left, right> pair, we use a 64 bit int whose first 32 bits are left and last 32 bits are right.
// This is safe because OpenFST state IDs are 32 bit ints.
inline std::uint64_t ProductKey(FusedStateId left, FusedStateId right) {
    return (static_cast<std::uint64_t> (static_cast<std::uint32_t>(left)) << 32) | static_cast<std::uint32_t>(right);
}

inline FusedStateId LeftState(std::uint64_t key) {
    return static_cast<FusedStateId>(key >> 32);
}

inline FusedStateId RightState(std::uint64_t key) {
    return static_cast<FusedStateId>(key & 0xffffffffULL);
}


inline bool IsFinite(const FusedWeight& weight) {
    // FusedWeight::Zero() in tropical semiring means +inf
    return weight != FusedWeight::Zero() && weight.Member();
}


// float32 loses precision once arc counts get large enough that accumulated
// weights lose bits. We hit this concretely: length-encoded weights of 0.1
// (float32) summed over many arcs, then pruned against an exact 0.0 threshold,
// dropped every final state because accumulated rounding error pushed the sum
// off the boundary. Rather than introduce a new OpenFST weight type, we keep
// all FST weights as float32 and only widen to double for this internal
// threshold/accumulation math.
constexpr double kDZero = std::numeric_limits<double>::infinity();
constexpr double kDOne = 0.0;

inline double ToDouble(const FusedWeight& weight){
    return IsFinite(weight)? static_cast<double>(weight.Value()) : kDZero;
}

inline FusedWeight ToWeight(double value){
    return value == kDZero? FusedWeight::Zero() : FusedWeight(static_cast<float>(value));
}

// Redefine tropical semiring operations in double
inline bool DIsFinite(double value) {return value != kDZero;}
inline double DPlus(double a, double b) {return std::min(a, b);}
inline double DTimes(double a, double b) {
    return (a==kDZero || b==kDZero)? kDZero : a + b;
}
inline bool DLess(double a, double b) {return a < b;}

// Exact distance-to-final calculation for one topologically sorted operand.
// These two small vectors form an admissible and consistent A* heuristic for
// the product: independently finishing the left and right paths can never cost
// more than finishing them while also requiring identical label sequences.
inline std::vector<double> AcyclicDistanceToFinal(const FusedFst& input) {
    const FusedStateId num_states = fst::CountStates(input);
    std::vector<double> distance(num_states, kDZero);
    for (FusedStateId state = num_states; state-- > 0;) {
        double best = ToDouble(input.Final(state));
        for (fst::ArcIterator<FusedFst> arcs(input, state); !arcs.Done(); arcs.Next()) {
            const FusedArc& arc = arcs.Value();
            const double candidate = DTimes(ToDouble(arc.weight), distance[arc.nextstate]);
            best = DPlus(best,candidate);
        }
        distance[state] = best;
    }
    return distance;
}

// Calls fn(left_arc, right_arc) once for every ordinary intersection arc.
// Both operands are epsilon-free acceptors. The right matcher uses its sorted
// input labels to avoid testing the full Cartesian product of outgoing arcs.
template <class Function>
inline void ForEachMatchingArc(const FusedFst& left, FusedStateId left_state,
                               fst::SortedMatcher<FusedFst>* right_matcher,
                               FusedStateId right_state, Function&& fn) {
    right_matcher->SetState(right_state);
    for (fst::ArcIterator<FusedFst> left_arcs(left, left_state); !left_arcs.Done(); left_arcs.Next()) {
        const FusedArc& left_arc = left_arcs.Value();
        if (!right_matcher->Find(left_arc.ilabel)) continue;
        for (;!right_matcher->Done(); right_matcher->Next()) {
            fn(left_arc, right_matcher->Value());
        }
    }
}

struct ProductNode {
    std::uint64_t key;
    double forward = kDZero;
    double backward = kDZero;
    double heuristic = kDZero;
    bool settled = false;
};

struct QueueItem {
    double estimate;
    double forward;
    std::uint64_t node;
};

struct QueueItemGreater {
    bool operator() (const QueueItem& left, const QueueItem& right) const {
        if (left.estimate != right.estimate) {
            return left.estimate > right.estimate;
        }
        if (left.forward != right.forward) {
            return left.forward > right.forward;
        }
        return left.node > right.node;
    }
};

inline void ValidateFusedOperand(const FusedFst& input, const std::string& name) {
    constexpr std::uint64_t required = fst::kAcceptor | fst::kAcyclic | fst::kTopSorted | fst::kILabelSorted | fst::kNoIEpsilons;
    const std::uint64_t properties = input.Properties(required, true);
    if ((properties & required) != required) {
        throw std::invalid_argument(
            name + " must be an epsilon-free, input-label-sorted, acyclic, "
                   "topologically sorted acceptor");
    }
}
} // namespace internal


// Intersects two MEDICC upper-pass acceptors and directly materializes only 
// paths within weight_threshold of the global optimum. Unlike
// fst::Intersect followed by prune_acyclic_std, this never stores the enormous
// raw product's arcs. Product states are discovered with A*, exact backward 
// product distances are then computed by regenerating matching arcs, and a 
// final traversal writes the small retained graph.
inline fst::script::MutableFstClass* intersect_prune_acyclic_std(
    const fst::script::FstClass& left_class,
    const fst::script::FstClass& right_class,
    float weight_threshold,
    float /*delta*/) {
    
    using namespace internal;

    if (left_class.ArcType() != FusedArc::Type() || 
        right_class.ArcType() != FusedArc::Type()) {
            throw std::invalid_argument(
                "intersect_prune_acyclic requires OpenFST standard-arc inputs (tropical semiring).");
    }

    const FusedFst* left = left_class.GetFst<FusedArc>();
    const FusedFst* right = right_class.GetFst<FusedArc>();
    if (left == nullptr || right == nullptr) {
        throw std::runtime_error("could not obtain typed OpenFST operands.");
    }

    ValidateFusedOperand(*left, "left operand");
    ValidateFusedOperand(*right, "right operand");

    if (!std::isfinite(weight_threshold) || weight_threshold < 0.0f) {
        throw std::invalid_argument("weight_threshold must be finite and nonnegative.");
    }
    const double threshold = static_cast<double>(weight_threshold);

    auto output = std::make_unique<fst::VectorFst<FusedArc>>();
    output->SetInputSymbols(left->InputSymbols());
    output->SetOutputSymbols(left->OutputSymbols());
    if (left->Start() == fst::kNoStateId || right->Start() == fst::kNoStateId) {
        return new fst::script::VectorFstClass(std::move(output));
    } // if either one is empty just return an empty FST 

    const auto left_distance = AcyclicDistanceToFinal(*left);
    const auto right_distance = AcyclicDistanceToFinal(*right);
    const auto heuristic_for = [&](FusedStateId left_state, 
                                   FusedStateId right_state) {
        return DTimes(left_distance[left_state], right_distance[right_state]);
    };

    const double start_heuristic = heuristic_for(left->Start(), right->Start());
    if (!DIsFinite(start_heuristic)) {
        return new fst::script::VectorFstClass(std::move(output));
    }

    std::vector<ProductNode> nodes; 
    nodes.reserve(1'000'000);
    std::unordered_map<std::uint64_t, std::uint64_t> node_by_key;
    node_by_key.reserve(1'000'000);
    std::priority_queue<QueueItem, std::vector<QueueItem>, QueueItemGreater> queue; 

    const std::uint64_t start_key = ProductKey(left->Start(), right->Start());
    nodes.push_back(ProductNode{start_key, kDOne, kDZero, start_heuristic, false});
    node_by_key.emplace(start_key, 0);
    queue.push(QueueItem{start_heuristic, kDOne, 0});

    double best_path = kDZero;
    double limit = kDZero;
    fst::SortedMatcher<FusedFst> right_matcher(right, fst::MATCH_INPUT);

    while (!queue.empty()) {
        const QueueItem item = queue.top(); 
        queue.pop();
        ProductNode& node = nodes[item.node];
        if (node.settled || item.forward != node.forward) continue;

        const double estimate = DTimes(node.forward, node.heuristic); 
        if (DIsFinite(best_path) && DLess(limit, estimate)) break;
        node.settled = true;

        // Keep a value copy across the callback, Discovering a new product 
        // state may grow `nodes` and invalidate the `node` reference. 
        const double current_forward = node.forward;

        const FusedStateId left_state = LeftState(node.key);
        const FusedStateId right_state = RightState(node.key);
        const double final_weight = DTimes(ToDouble(left->Final(left_state)), ToDouble(right->Final(right_state)));
        const double final_candidate = DTimes(current_forward, final_weight);
        if (DLess(final_candidate, best_path)) {
            best_path = final_candidate;
            limit = DTimes(best_path, threshold);
        }

        ForEachMatchingArc(*left, left_state, &right_matcher, right_state,
            [&](const FusedArc& left_arc, const FusedArc& right_arc) {
                const double arc_weight =
                    DTimes(ToDouble(left_arc.weight), ToDouble(right_arc.weight));
                const double candidate = DTimes(current_forward, arc_weight);
                const double heuristic = heuristic_for(left_arc.nextstate, right_arc.nextstate);
                if (!DIsFinite(heuristic)) return;
                if (DIsFinite(best_path) && DLess(limit, DTimes(candidate, heuristic))) return; // based on heuristic this node will not be inside the pruned graph.

                const std::uint64_t key = ProductKey(left_arc.nextstate, right_arc.nextstate);
                auto found = node_by_key.find(key);
                std::uint64_t next;
                if (found == node_by_key.end()) {
                    // new node not in the nodes visited, add this node to nodes
                    if (nodes.size() >= std::numeric_limits<std::uint64_t>::max()) {
                        throw std::overflow_error("too many fused product states");
                    }
                    next = static_cast<std::uint64_t>(nodes.size());
                    nodes.push_back(ProductNode{key, candidate, kDZero, heuristic, false});
                    node_by_key.emplace(key, next);
                    queue.push(QueueItem{DTimes(candidate, heuristic), candidate, next});
                } else {
                    // node already in the nodes visited, compare and see if it needs updating
                    next = found->second;
                    ProductNode& next_node = nodes[next];
                    if (DLess(candidate, next_node.forward)) {
                        if (next_node.settled) {
                            throw std::runtime_error("Inconsistent product heuristic settled a state too early.");
                        }
                        next_node.forward = candidate;
                        queue.push(QueueItem{DTimes(candidate, next_node.heuristic), candidate, next});
                    }
                }
            });
    }

    if (!DIsFinite(best_path)) {
        return new fst::script::VectorFstClass(std::move(output));
    }

    // Only settled states whose independent lower bound is within the final
    // threshold can possibly occur on a retained product path.
    std::vector<std::uint64_t> active;
    active.reserve(nodes.size());
    for (std::uint64_t index = 0; index < nodes.size(); ++index) {
        const ProductNode& node = nodes[index];
        if (node.settled && !DLess(limit, DTimes(node.forward, node.heuristic))) {
            active.push_back(index);
        }
    }
    std::sort(active.begin(), active.end(), [&](std::uint64_t a, std::uint64_t b) {
        return nodes[a].key > nodes[b].key;
    });

    // Arc destinations have larger component state IDs, hence larger packed
    // keys. Descending key order is therefore a valid reverse topological order
    // for exact product distance to final dynamic programming
    for (const std::uint64_t index : active) {
        ProductNode& node = nodes[index];
        const FusedStateId left_state = LeftState(node.key);
        const FusedStateId right_state = RightState(node.key);
        double backward = DTimes(ToDouble(left->Final(left_state)), ToDouble(right->Final(right_state)));
        ForEachMatchingArc(*left, left_state, &right_matcher, right_state,
            [&](const FusedArc& left_arc, const FusedArc& right_arc) {
               const auto found = node_by_key.find(ProductKey(left_arc.nextstate, right_arc.nextstate));
                if (found == node_by_key.end()) return;
                const ProductNode& next_node = nodes[found->second];
                if (!next_node.settled || !DIsFinite(next_node.backward)) return;
                const double candidate = DTimes(DTimes(ToDouble(left_arc.weight), ToDouble(right_arc.weight)),
                                                next_node.backward);
                backward = DPlus(backward, candidate);
            });
        node.backward = backward;
    }

    const auto start_found = node_by_key.find(start_key);
    if (start_found == node_by_key.end() || !DIsFinite(nodes[start_found->second].backward)) {
        return new fst::script::VectorFstClass(std::move(output));
    }

    std::vector<FusedStateId> output_state(nodes.size(), fst::kNoStateId);
    std::queue<std::uint64_t> pending;
    const std::uint64_t start_index = start_found->second;
    output_state[start_index] = output->AddState();
    output->SetStart(output_state[start_index]);
    pending.push(start_index);

    while (!pending.empty()) {
        const std::uint64_t index = pending.front();
        pending.pop();
        const ProductNode& node = nodes[index];
        const FusedStateId left_state = LeftState(node.key);
        const FusedStateId right_state = RightState(node.key);

        const double final_weight = DTimes(ToDouble(left->Final(left_state)), ToDouble(right->Final(right_state)));
        if (!DLess(limit, DTimes(node.forward, final_weight))) {
            output->SetFinal(output_state[index], ToWeight(final_weight));
        }

        ForEachMatchingArc(*left, left_state, &right_matcher, right_state,
            [&](const FusedArc& left_arc, const FusedArc& right_arc) {
                const auto found = node_by_key.find(ProductKey(left_arc.nextstate, right_arc.nextstate));
                if (found == node_by_key.end()) return;
                const std::uint64_t next = found->second;
                const ProductNode& next_node = nodes[next];
                if (!DIsFinite(next_node.backward)) return;
                const double arc_weight = DTimes(ToDouble(left_arc.weight), ToDouble(right_arc.weight));
                const double through_arc = DTimes(DTimes(node.forward, arc_weight), next_node.backward);
                if (DLess(limit, through_arc)) return;

                if (output_state[next] == fst::kNoStateId) {
                    output_state[next] = output->AddState();
                    pending.push(next);
                }
                output->AddArc(output_state[index], FusedArc(left_arc.ilabel, left_arc.ilabel, ToWeight(arc_weight), output_state[next]));
            });
    }

    fst::Connect(output.get());
    return new fst::script::VectorFstClass(std::move(output));

}
    
} // namespace medicc_fstlib_extension

#endif // MEDICC_FSTLIB_EXTENSION_ACYCLIC_PRUNE_H_