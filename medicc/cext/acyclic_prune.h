#ifndef MEDICC_FSTLIB_EXTENSION_ACYCLIC_PRUNE_H_
#define MEDICC_FSTLIB_EXTENSION_ACYCLIC_PRUNE_H_

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <fst/expanded-fst.h>
#include <fst/fstlib.h>
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

} // namespace medicc_fstlib_extension

#endif // MEDICC_FSTLIB_EXTENSION_ACYCLIC_PRUNE_H_