Exact fallback for expensive ancestor joins
==========================================

A few capped-model upward joins explore tens or hundreds of millions of history
states although their final ancestor languages are small. The default join now
reuses the existing chromosome-factorization prototype only after a fixed-WGD
search expands 250,000 states. This threshold selects an algorithm; it never
limits the answer, drops a state for being expensive, or changes prune-zero ties.

Read ``medicc/full_fst/model.py`` for scalar witness bounds, then
``medicc/cext/full_fst/ancestors/src/nfa_engine.hpp`` for the trigger and
``segmented_join.hpp`` for the exact fallback. The latter adapts experiment 24's
existing implementation; it uses the same Boost.Graph search and OpenFST language
operations as the ordinary join.

Safe cuts and costs
-------------------

A cut is allowed only immediately after an explicit X where both child DFAs
have one state. Every primitive must reset to its start state at X, and that
start must have zero final cost. These checks preserve the child languages'
correlations. If no safe cut exists, the original search continues unchanged.

At a safe cut, prefix and suffix choices are independent. For one fixed pair of
WGD counts, the optimum is the sum of the independent segment optima plus the
child offsets and one whole-genome WGD charge. Every globally optimal word must
use an optimal word in every segment; concatenating all segment-optimal languages
therefore preserves the complete optimal language, including ties.

The restart releases the original explored graph before allocating segment
searches. Small segments are solved first to tighten the remaining budget.
Admissible bounds reserve the cost of unsolved segments. Answers are stored by
original position and concatenated in chromosome order, so evaluation order does
not reorder the genome. A segment never chooses its own WGD count: every segment
uses the pair selected by the outer join. The local search's WGD charge is removed
from each segment result and charged once in the final whole-genome offset.

Segment searches retain the existing chromosome and sibling bounds. They cannot
restart recursively. The outer join still keeps every fixed-WGD branch tying at
the best score. No new two-WGD sibling abstraction is introduced in this change.

Other reused optimizations
--------------------------

The backward pass uses the already computed full forward heuristic when deciding
which expanded states can participate in an optimum. It still retains equality.
This avoids regenerating transitions from states ruled out by the tighter bound.

``Model.bound_scorer`` uses ``direct_pair`` for the exact scalar witness bound,
even when ``Model.scorer`` uses the independent reference backend. The reference
backend remains selected for final scores when requested by the workflow. Bounds
do not select ancestors; the native join/resolve still makes the same choices.
When both scorers are direct they share one object. The model closes an additional
bound scorer separately when it owns one.

Validation and limitations
--------------------------

Experiment 27, ``27-capped-ancestor-outliers``, contains the orchestration scripts.
Data, frozen child messages, logs and results remain outside this tool folder in
``results_treebench_development/medicc_capped_ancestor_opt_20261016``.

The normal build uses the 250,000-state trigger. A separate test-only build with
``-DMEDICC_ANCESTOR_SEGMENT_TRIGGER=1`` forces the fallback on small cases. The
weighted oracle and boundary tests cover factorized and correlated languages,
partial cuts, offsets, final weights, mixed separator positions, leading/trailing
and consecutive separators, plus serial/parallel WGD branches. The test-only
library is never installed into the production checkout.

``tests/check_ancestor_segments.py TOOL OUTPUT [--expected BASELINE] [--forced]``
compares complete weighted languages and checks fallback coverage in the forced
build. Other existing capped, weighted-join, binding, tie and API checks remain
applicable. Complete-tree comparisons must preserve profiles, branch scores,
distance matrices and accuracy metrics before a release is used for production.

This is an exact optimization, not a general runtime guarantee. Correlated child
languages may have no safe cut, and an individual segment can still be expensive.
The trigger also incurs some repeated work before switching algorithms. Timing
claims should distinguish controlled fixture comparisons from historical runs on
other allocations.

Release evidence
----------------

Controlled ABBA comparisons on three frozen weighted-language joins preserve
the complete languages and improve mean runtime by 1.13x, 2.09x and 2.55x.
The 500-leaf meso iteration11 truth reconstruction changes from a historical
1413.47 s to 115.32 s, with identical complete profiles, branch scores and Hamming
results. On the saved FastME topology, reconstruction changes from 1345.06 s to
163.86 s with identical profiles and branch scores. Peak process memory changes
from about 30 GiB to 5.5--5.7 GiB. These tree timings use different allocations,
so the controlled fixtures provide the stronger timing evidence.

The direct and reference pair libraries are byte-identical to the prior release.
Small FastME/truth workflows also reproduce their scientific outputs. The large
FastME test reuses its saved topology rather than recomputing unchanged distances.
