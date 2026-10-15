A smaller sibling problem guides the exact search
================================================

The expensive join must find an ancestor profile that can produce **both** child
profiles. The previous suffix bound treated the children independently. It could
therefore regard two incompatible explanations as cheap, leaving millions of
states for the exact search to examine.

The new bound solves a small, relaxed sibling problem first. It still requires
both sides to consume the **same ancestor symbol at each column**. Its answer is
a lower bound used by the existing Boost A* search. The full model, event costs,
and set of optimal ancestor profiles are unchanged.

The small model
---------------

``A`` is the fused loss/LOH/gain transducer; ``D`` is capped doubling. The exact
one-WGD branch is represented as ``AD, A``, using the OpenFST-reduced ``AD`` factor.
For the bound only, each primitive state is placed in one of four classes:

* neither a free decrease nor a free increase is available;
* a free decrease is available;
* a free increase is available;
* both are available.

Here "free" means that an outgoing arc has weight zero. For ``A``, changes are
measured relative to the input copy number. For ``AD``, they are measured relative
to ``min(2 * input, 8)``. Every original arc is projected to its two state classes;
duplicate projected arcs keep their minimum weight. The classification controls
how informative the bound is. Its safety comes from projecting **all** arcs.

A zero-WGD side has four abstract states. A one-WGD side has two four-state
factors, giving sixteen states. The paired lookup table therefore has:

.. list-table:: Abstract states per column
   :header-rows: 1

   * - Left/right WGD counts
     - Paired states
   * - 0 / 0
     - 16
   * - 0 / 1 or 1 / 0
     - 64
   * - 1 / 1
     - 256

Branches involving **two WGDs** continue to use the previous certified bound and
exact search. They are not discarded or approximated. The constructor rejects
unsupported counts if called incorrectly.

The computation, in Python-like notation
---------------------------------------

For each column, the child's mask allows every symbol appearing in that column
of its retained language. This deliberately ignores correlations between columns
and nonnegative child-message weights.

.. code-block:: python

   # s and t each represent a pair of abstract primitive states.
   remaining[length, :] = 0
   for column in reversed(range(length)):
       for s in paired_states:
           remaining[column, s] = min(
               left_cost(s.left, t.left, ancestor_symbol, left_mask[column])
               + right_cost(s.right, t.right, ancestor_symbol, right_mask[column])
               + remaining[column + 1, t]
               for ancestor_symbol in alphabet
               for t in paired_states
           )

   lower_bound = max(previous_bound, remaining[column, abstract_state])
   if prefix_cost + lower_bound > incumbent:
       discard_this_search_state()

``left_cost`` and ``right_cost`` minimize over the allowed child symbols. The
shared ``ancestor_symbol`` is what couples the two explanations. Constant child
offsets and WGD charges remain accounted for by the existing search, outside this
table. Transition matrices are cached by the two column masks; at most 64 joint
matrices are retained. Eviction only causes recomputation.

Why this cannot remove an optimum or an optimal tie
--------------------------------------------------

Every real suffix maps to a path through the abstract classes with the same
ancestor and child labels. Its projected arc weights are no larger. Its child
labels remain allowed by the column masks, and dropping nonnegative child/final
weights can only lower its cost. The abstract optimum is therefore no greater
than the true remaining cost.

Taking the **maximum** with the previous bound strengthens it without counting
any cost twice. A state is discarded only when its best possible total is
**strictly greater** than the feasible incumbent. Equality remains eligible, so
the existing backward pass still collects every optimal ancestral word.

Reading the C++
---------------

The vectors are ordinary owned arrays. Expressions such as ``q * 4 + r`` pack two
small state numbers into one array index. ``side_state`` decodes the existing
native tuple and maps it to these class numbers. The two pointers named ``next``
and ``current`` merely view adjacent rows of the owned ``values`` vector; they do
not allocate or own memory. ``const ... &`` borrows a table without copying it.

This follows the established *pattern database* method: compute exact distances
in an abstract problem and use them as admissible search heuristics. The
MEDICC-specific adaptation is the weighted FST quotient and the shared-ancestor
constraint. See Culberson and Schaeffer, `Pattern Databases (1998)
<https://webdocs.cs.ualberta.ca/~jonathan/publications/ai_publications/compi.pdf>`_,
especially its lower-bound lemma and maximum combination of bounds.

Exact OpenFST factor reduction
------------------------------

The production search now represents A (D A)^k as (AD)^k A. The generated AD
factor is the exact weighted composition of A and capped D, reduced with
OpenFST paired-label encoding, determinization and minimization. It has 74 states
and 3,108 arcs, compared with A's 212 states and 15,688 arcs. We do not build the
full two-WGD history relation. ``scripts/generate_capped_ad.py`` reproduces the
asset, verifies weighted equivalence after minimization and canonical text export,
and checks zero-cost chromosome resets. No new runtime dependency is introduced.

This compact factor is used for bottom-up joins and pair distances. Downward
reconstruction deliberately retains the original A,D,A factor order, because
its historical choice among equal-cost profiles depends on primitive state
ordering. Two nine-tree trials exposed that difference: full retained DFAs were
byte-identical, with the same conditional costs, but a different optimal word
was chosen. The final implementation preserves the previous representative.
``tests/check_resolve_tie.py`` protects a small counterexample; an independent
5,000-case weighted-DFA audit matched both costs and chosen words.

Suffix-table preparation
------------------------

The existing exact last-A suffix recurrence is unchanged. Its arcs are indexed
by the observed child output label before the dynamic program, so a sparse
child message visits only compatible arcs. Input parity restrictions, child
weights, final weights and explicit X resets are preserved. This removes work
from setup after the search explosion is reduced.

Measured full-tree results
--------------------------

Nine selected difficult 250-leaf trees used eight workers. The first completed
integrated trial for each tree is listed below. Previous capped timings use
commit 9fe5107 (including the earlier warm starts); strict timings are the
preserved uncapped pilot. These comparisons use separate cluster allocations.
The strict implementation has not received these new optimizations, so this is
an implementation comparison, not a claim that the capped model is inherently
cheaper. All final profiles and all branch costs are identical to the preserved
capped results, including the two cases that exposed the traceback tie issue.

.. list-table:: Ancestral reconstruction, seconds
   :header-rows: 1

   * - Dataset / iteration
     - Previous capped
     - New capped
     - Old uncapped
   * - meso_low_wgd / 12
     - 105.02
     - 26.45
     - 40.63
   * - meso_low_wgd / 7
     - 559.64
     - 25.36
     - 35.84
   * - meso_low_wgd / 9
     - 107.05
     - 27.74
     - 39.61
   * - meso_low_wgd / 17
     - 200.94
     - 33.14
     - 33.78
   * - ucec_mid_wgd / 8
     - 1162.39
     - 55.29
     - 47.67
   * - pancancer / 18
     - 95.08
     - 25.65
     - 35.39
   * - pancancer / 9
     - 108.68
     - 20.01
     - 37.32
   * - pancancer / 2
     - 218.14
     - 33.77
     - 41.55
   * - lgg_high_wgd / 9
     - 154.43
     - 20.29
     - 33.85

Median per-tree speedup over the previous capped implementation is 6.06x;
median new reconstruction time is 26.45 seconds. The median new/uncapped ratio
is 0.708, with range 0.536--1.160. The hardest tree now takes 55.29 seconds,
compared with 1,162.39 seconds before and 47.67 seconds uncapped. Its peak RSS is
about 2.40 GiB, close to the old uncapped run and below the earlier capped 9.5 GiB.

A separate same-allocation comparison against the immediate parent commit
eb96a65 changed meso iteration 12 from 107.59 to 25.29 seconds (4.25x), with
identical outputs. Frozen heavy joins also preserve their complete weighted
ancestor-language signatures: meso iteration 7 changed 529.90 to 7.09 seconds;
ucec iteration 8 changed 1,082.36 to 10.52 seconds. Total explored states fell
from 24.28 million to 293,498 and from 58.53 million to 366,659 respectively.

Pairwise distances remain the main cost
--------------------------------------

Final integrated source was also checked on two complete 250-leaf distance
matrices (31,375 pairs each), using eight workers. Every distance is identical
to the previous capped matrix. These are separate allocations, not paired
same-node timing trials.

.. list-table:: Distance calculation, seconds
   :header-rows: 1

   * - Dataset / iteration
     - Previous capped
     - New capped
     - Old uncapped
     - New / uncapped
   * - meso_low_wgd / 9
     - 261.23
     - 134.49
     - 83.22
     - 1.62x
   * - ucec_mid_wgd / 8
     - 4007.93
     - 577.21
     - 129.48
     - 4.46x

Peak RSS fell from 1,359 to 892 MiB and from 7,238 to 3,623 MiB respectively.
The final source is 1.94x and 6.94x faster than the previous capped distance
implementation, but it still has a substantial cost relative to uncapped on
UCEC. Earlier intermediate-factor timings of 112.81 and 487.36 seconds are
not the final integrated measurements.

FastME RF distances remain 2 for both capped matrices, versus 10 and 4 for
the respective uncapped runs. This preserves the pilot's apparent accuracy
benefit; two selected cases do not establish a general accuracy improvement.

A paired ABBA trial of the sibling pattern bound in the distance engine
was 21.8% slower and was rejected. Profiling one remaining expensive pair
placed 3.04 of 3.66 seconds in the two-WGD/one-WGD branch. A final bounded
OpenFST attempt to reduce AD composed with AD produced 4,603 states and
782,510 arcs before determinization, but determinization exceeded 12 GiB
at 3 minutes 2 seconds. No such factor is installed. The two-WGD branch is
therefore the main remaining optimization target, not a solved bottleneck.

Validation and scope
--------------------

The release candidate passed the weighted all-tie oracle, literal FST oracle,
high-CN tests, bindings, 64 additional weighted chromosome-language joins,
18 boundary/language cases, 5,000 conditional traceback comparisons, and all
nine full-tree output comparisons. The three frozen real joins compare their
entire weighted DFA languages, not just one selected ancestor. The AD generator
also verifies its weighted export using OpenFST.

Exact segmentation at singleton chromosome boundaries, whole-branch witness
prepasses, and a cost-first segmented search were tested separately. They were
correct but less effective on the hardest case and are not part of this change.
An adaptive probe/restart variant passed tests but was not selected: profiling
showed suffix setup was the more valuable next target.

These are selected hard cases, not an unbiased model-accuracy study. Two-WGD
branches retain their previous bounds; this does not establish a polynomial
worst-case search algorithm. The model, local prune_weight=0 rule and whole-genome
WGD accounting are unchanged. RF and ancestral accuracy cannot change when all
distances, profiles and branch costs are preserved.

Reproduction
------------

Source-only scripts are in
``pipeline-simulation-benchmarking/chenxi_experiments/24-capped-search-explosion``.
Inputs, immutable snapshots, complete messages, timing records and logs remain
outside this checkout under
``/home/cnie/scratch_folder/results_treebench_development/medicc_search_explosion_20261015``.
All task jobs use ``preempt-schwarzr``. ``collect.py`` reports completed frozen
joins; ``collect_trees.py`` reports completed full-tree output checks.
