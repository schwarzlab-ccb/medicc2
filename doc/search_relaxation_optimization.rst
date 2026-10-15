Avoiding redundant forward relaxations
======================================

The additional change
---------------------

After the feasible-path warm start, the largest joins still generate many
histories reaching states whose prefix cost is already as good or better.
Boost.Graph correctly rejects these relaxations. The MEDICC adapter now checks
that same condition before constructing and handing the edge to Boost::

    candidate_cost = source.prefix_cost + transition.cost
    if candidate_cost < destination.prefix_cost:
        forward_edges.append(transition)
    else:
        relaxation_skips += 1

This is in ``ancestors/src/nfa_engine.hpp``, in the callback passed to
``library_search::run``. Boost.Graph still owns A*, its queue and relaxation.
The model, pair-distance engine and top-down reconstruction are unchanged.
No additional library or cache is introduced.

Why equal-cost ancestors survive
-------------------------------

The forward search computes best prefix distances. It does not store the
ancestor language or its predecessor edges. Equal-cost routes to the same state
cannot improve that distance, so Boost would reject them too. Strictly better
routes still reach Boost, including routes reopening a previously visited state.

The backward pass is unchanged: it regenerates transitions from retained states,
computes exact suffix costs, and retains every tight labeled transition. Thus
multiple ancestor symbols reaching the same state with equal cost still survive.
This guard must not be copied into that backward language-building pass.
All cutoffs on the optimal score still retain equality.

The trace field ``relaxation_skips`` counts avoided forward relaxations.
``generated`` also counts backward and output work, so it is not a forward-only
denominator. On the difficult meso iteration 7 join, the filter avoids about
1.28 billion forward relaxations but leaves the 24,010,913 explored states
unchanged. Transition generation and backward work remain expensive.

Measured effect
---------------

Baseline is commit ``9fe5107`` including both previous warm-start optimizations.
Two sequential, same-allocation comparisons on the selected 250-leaf meso
iteration 12 tree used eight workers:

* Baseline then candidate: 97.98 -> 90.54 seconds (7.6% less time).
* Candidate then baseline: 96.55 -> 87.37 seconds (9.5% less time).

The reverse-order check guards against attributing every warm-cache advantage
to the candidate. These are two controlled trials on one selected tree, not a
universal speedup. Complete profiles and every branch cost matched in both.
A separate-allocation meso iteration 7 check changed 559.64 -> 539.01 seconds;
The largest ucec iteration 8 check changed 1162.39 -> 1145.35 seconds (about
1.5% less time). All profiles and branch costs matched on both trees.
These smaller gains should be interpreted cautiously because allocations differ.
This filter does not reduce the number of search states and promises no RAM
reduction. Whole-process RSS measurements are retained in the external report.

Validation and rejected experiments
-----------------------------------

The candidate passed 192 weighted all-tie joins against the literal oracle,
240 oracle messages, 240 resolutions, 560 directed edges, 64 incremental NNI
neighbors, serial/parallel NNI, capped high-CN regressions and binding checks.
Complete benchmark profiles and branch costs are compared with preserved outputs.
These language checks are stronger than comparing one selected ancestor alone.

Prepared-state and edge caches, one-column lookahead and a chromosome-suffix
bound were also tested. Their overhead was not compensated by useful reductions
on the pilot, so none is included. Moving this guard before the heuristic write
and removing its diagnostic counter also regressed in its controlled trial.

An OpenFST A D A fragment composed successfully to 29,585 states and 17,336,810
arcs in 1.03 seconds, but its determinization exceeded the 12 GiB trial limit.
This is a resource-limited result for that configuration, not a claim that
OpenFST cannot express the relation. The existing lazy factorization is retained.

Reproduction and model accuracy
-------------------------------

Source-only trial scripts live in
``pipeline-simulation-benchmarking/chenxi_experiments/23-capped-search-reuse``.
Snapshots, logs, input references and reports are outside the tool checkout at
``/home/cnie/scratch_folder/results_treebench_development/medicc_capped_reuse_20261015``.
Run that experiment's ``collect.py`` to collect completed results. All its jobs
use ``preempt-schwarzr``; pending runs are excluded from completed comparisons.

The unchanged pair engine also tested another selected 250-leaf meso tree
(iteration 9): strict doubling gave RF 10, capped doubling RF 2. Distance time
was 83.22 versus 261.23 seconds. Together with the previous iteration 12 result
(RF 14 -> 2), this motivates further evaluation of the capped model, but these
selected hard cases do not establish its accuracy on an unbiased dataset.
The search filter itself preserves scores and cannot improve RF or ancestors.

The additional ucec iteration 8 RF comparison was still running when this change
was committed (job ``6828463_4``). It is not included as a completed RF result
here. Collector job ``6829181`` refreshes the external report after it finishes.
