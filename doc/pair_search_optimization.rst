Capped-WGD pair-search optimization
==================================

This branch preserves D(c) = min(2*c, 8), full L_LOH, and at most two WGDs
per directed edge. It speeds up shared-ancestor pair distances. The Python
API, directed scorer, ancestral reconstruction, and all-tie message semantics
are unchanged.

What was expensive
------------------

The scorer fixes the left/right WGD counts, then solves chromosomes independently.
Previously, each chromosome received whatever budget remained under the whole
pair's incumbent score. This can be much larger than a feasible chromosome cost.
A* then allocates many expensive successor states before it reaches its first
complete path. In one profiled search, 24,579 expansions generated about
831 million transitions and retained almost 42 million states.

What changed
------------

``pairwise.cpp`` now builds feasible upper bounds in three steps:

1. ``chromosome_witness`` tries constant positive ancestors and a binwise
   ancestor. It applies the fixed mandatory doublings, then valid interval
   losses/gains. Each candidate is a realizable history.
2. For each WGD-count pair, the chromosome witnesses concatenate at X.
   Their costs sum, and the WGD charges are added once per genome. This supplies
   a better whole-pair incumbent without selecting different WGD counts for
   different chromosomes.
3. ``greedy_witness`` walks the existing primitive transitions for one
   chromosome. At each bin it tries a common ancestor symbol and chooses the
   best compatible edge on each child side, using the existing suffix heuristic.
   It avoids enumerating the two edge lists' Cartesian product. A complete
   walk gives another feasible upper bound. An unsuccessful walk keeps the
   previous bound.

The greedy walk is only a way to find a valid history cheaply. Its score is
never accepted as optimal without the exact search proving that nothing cheaper
exists. The logic is equivalent to::

    upper = cost_of_a_known_valid_history()
    upper = try_a_greedy_history_or_keep(upper)
    cutoff = min(remaining_genome_budget, upper)
    improvement = exact_search_for_cost_strictly_below(cutoff)
    if improvement.exists:
        return improvement.cost
    if upper < remaining_genome_budget:
        return upper  # feasible, and exact search proved no cheaper path exists
    return cannot_improve_this_genome_branch

The integer weights allow a further small optimization: score-only transition
generation uses ``remaining_budget - 1``, since a path tying the cutoff cannot
improve the score. Searches constructing ancestral languages still retain all
paths tying their optimum.

Implementation scope
--------------------

The runtime change is in two existing C++ files:

* ``medicc/cext/full_fst/direct_pair/src/pairwise.cpp``: feasible witnesses,
  chromosome upper bounds, and the whole-genome incumbent.
* ``medicc/cext/full_fst/direct_pair/src/nfa_engine.hpp``: early strict-cutoff
  pruning for score-only callers.

The code reuses the existing transition enumeration, suffix potentials,
state storage, and priority queue. It adds no runtime dependency, approximate
search cutoff, beam, or state-count limit. OpenFST's existing responsibilities
are unchanged.

What was tried and rejected
---------------------------

The experiments first tightened lower bounds. One retained boundary-relaxation
information through the chromosome suffix; another additionally required
neighboring boundaries to agree on their common ancestral copy number.
OpenFST composition and SciPy shortest paths generated safe credit tables for
already-open events in the compressed FST states. Both passed correctness
checks, but their extra work outweighed the pruning benefit on the measured row:
about 0.61x and 0.88x baseline speed respectively. Neither is included here.

A raw, uncompressed primitive composition was also slower and was stopped.
Caching lower bounds from exhausted searches yielded only about 1.05x in the
row trial and was not adopted. Early strict-cutoff pruning alone gave about
1.21x; simple local witnesses alone gave about 1.35x. The large improvement
came from the feasible-path warm start used with the exact search.

Validation and measurements (2026-10-15)
---------------------------------------

The optimized build passed:

* 52,488 exhaustive short-string comparisons across the two pair backends.
* 960 longer-profile comparisons, 12 targeted high-CN cases, 24 all-tie joins,
  and 18 conditional resolves.
* 14,745 literal fixed-WGD witness/cutoff checks in ``tests/check_pair_bounds.py``.
  This test derives each fixed-WGD relation independently from L_LOH, G and
  mandatory capped doubling, and compiles its test-only probe in a temporary
  directory. It tests witnesses and cutoffs below, at, and above the optimum.
* Public API and binding/ownership checks.
* Identical scores on a 249-pair real-data row and 54 pairs spread across all
  nine frozen difficult 250-leaf cases.

For the high-WGD LGG iteration-9 tree, with eight workers:

=============================== ============= =============
Measurement                     Original      Optimized
=============================== ============= =============
Distance-matrix time             4,162.73 s    79.08 s
Unique pairwise distances        31,375        31,375
Distance differences             --            0
Inferred topology file           baseline      identical
=============================== ============= =============

That is a 52.64x speedup for this case. The optimized distance run peaked at
about 536 MiB RSS. The second full matrix, meso iteration 12, completed in
329.19 seconds; its 249 previously measured row entries match exactly.
Its original capped full matrix had not completed at the time of this note.

These are single runs on a shared cluster, not controlled repeated timing
estimates. The nine-case spread contains six fixed pairs per case, not nine
complete matrix comparisons. This optimization preserves the capped model's
inferred results; it does not establish that capped WGD is biologically more
accurate than strict WGD, nor accelerate the separate ancestral-reconstruction
stage.

Reproduce
---------

Build normally, then run::

    python -B tests/check_pairs.py
    python -B tests/check_capped.py
    python -B tests/check_pair_bounds.py
    python -B tests/check_api.py
    python -B tests/check_bindings.py

The benchmark repository contains scripts in
``chenxi_experiments/21-capped-wgd-pilot``. ``collect_optimization.py`` collects
paired timings and score agreement; ``optimization_cases.sbatch`` runs the
nine-case spread; ``optimization_matrix.sbatch`` runs full distance matrices
and FastME. All use the task-specific ``preempt-schwarzr`` QoS. Experimental
inputs, outputs, candidate builds, profiling logs, and reports remain outside
this tool checkout.
