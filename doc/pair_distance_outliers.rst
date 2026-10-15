Restore short-chromosome-first distance evaluation
==================================================

This branch changes only the order in which the fixed-string pair scorer solves
chromosome factors. It restores an existing optimization from the original
``codex_effort/lloh_optimized/engine.cpp`` (``AR_ORDER``, around line 220), also
present in its lazy scorer. The original derivation is in
``codex_effort/reference_check/budgeted/README.md`` under "Chromosome order and
lazy outer A". These sources live under experiment 19's
``scripts/kernel_speedup/full_wgd`` directory.

The later factorized scorer, imported in f050793 and inherited by the capped
release fbc0a2d, instead sorted chromosomes by decreasing lower bound. The
ordering policy was therefore different even though chromosome factorization,
feasible incumbents and global WGD enumeration were already present. This is a
restoration of a known policy, not a new composition engine or a new model.

Why the outliers were expensive
------------------------------

The meso_low_wgd 250-leaf iteration19 diagnostic found shortened chromosomes of
only 35--40 columns that generated over 100 million candidate transitions in
one-WGD/two-WGD searches. Millions of product states accumulated even though
these WGD combinations could ultimately lose to another combination.

A large lower bound is not necessarily cheap to make exact. Solving that
chromosome first can spend substantial work before inexpensive chromosomes
have supplied enough evidence to reject the global WGD choice. Solving short
chromosomes first makes that evidence available earlier. Length here means the
pair-shortened chromosome length, not the original number of genomic bins.

The production progress logger counts pairs only after complete matrix rows
return. Its reported progress is a lower bound on the number of individual
pairs already scored inside workers. Long gaps between row completions cannot
be extrapolated as a reliable per-pair throughput estimate.

Why the answer is unchanged
---------------------------

For fixed WGD counts k_left and k_right, the scorer computes

.. math::

   k_\mathrm{left} + k_\mathrm{right}
      + \sum_c d_c(k_\mathrm{left}, k_\mathrm{right}).

Chromosome separators reset primitive event states. The chromosome costs are
independent summands for these fixed counts. The code still minimizes over the
same nine global WGD choices, charging WGDs once for the whole genome.

During evaluation, exact costs replace lower bounds. If their sum plus the
remaining lower bounds reaches an already feasible whole-genome score, this
WGD choice cannot improve the answer. Permuting summands changes the amount of
work before that proof is obtained, but cannot change its validity. No threshold,
beam, state cap, alphabet, event weight, or feasible history is changed.

Only scalar undirected pair scoring uses the new order. Directed scoring,
ancestral reconstruction, the independent reference scorer, and ancestor tie
selection are unchanged. Ancestral and reference binary hashes were checked
against fbc0a2d. There are no new dependencies.

Validation and results
----------------------

* Existing high-CN literal-FST regression: 960 pair comparisons, 24 complete
  ancestral-language joins, 18 conditional resolutions and 12 targeted checks.
* New regression: 36 independent reference cases with 216 checks of cold/warm
  caches, swapped children and permuted chromosomes, including empty separators.
* Existing Python/native ownership and binding regression passed.
* Three completed full matrices exactly match the production capped outputs:
  94,125 distinct pair distances in total.

.. list-table:: Complete distance matrices, eight workers
   :header-rows: 1

   * - 250-leaf dataset
     - Production capped seconds
     - Restored ordering seconds
     - Equal matrix
   * - meso iteration2
     - 3654.52
     - 110.68
     - Yes
   * - meso iteration9
     - 112.97
     - 40.99
     - Yes
   * - UCEC iteration8
     - 492.29
     - 64.54
     - Yes

These whole-matrix runs used different cluster allocations. They are observed
wall times, not a controlled hardware-only speedup claim. A separate repeated
baseline/candidate/candidate/baseline trial uses one allocation and reports
native scoring work on 12 frozen outlier pairs and 54 earlier control pairs.
The selected outlier pairs took 142.50 / 142.68 seconds in the two baseline
passes and 0.818 / 0.828 seconds in the restored-order passes, approximately
173x faster. Expanded states fell from 9,632,959 to 204,685 per pass. The earlier
54-pair spread took 3.052 / 3.069 seconds versus 0.679 / 0.681 seconds, about
4.5x faster. All scores match in every pass. These are sums of native scoring
times, excluding input loading, and are not full-workflow speedups. The 12
outlier pairs are a selected development set.

The full iteration19 matrix completed in 420.31 seconds with 2,247 MiB peak RSS.
Its production counterpart was still running after over 77 minutes, so a full
production-matrix comparison is pending; do not claim all 31,375 entries have
been independently matched. The sampled iteration19 outlier pairs do match. Audit job 6839050 is queued
with an after-success dependency on the original production job; it will
write ``full_matrix_audit.json`` after comparing all four complete matrices.
The former uncapped distance time was 862.04 seconds on a different allocation.

For iteration2, process peak RSS fell from 12,312 MiB in production to 763 MiB
in this trial. These observations demonstrate a workload-specific improvement;
no universal upper bound on capped-model runtime or memory is claimed.

The running 900-dataset experiment uses its frozen original snapshot. It was
not restarted or modified. Its existing results and timings remain valid for
that snapshot; the new timings must be labeled separately.

Reading and reproduction
------------------------

The implementation change is the ``sort`` comparator in
``medicc/cext/full_fst/direct_pair/src/pairwise.cpp``. Read its four-line comment,
then the surrounding global-WGD loop. The regression is
``tests/check_pair_chromosomes.py``.

Scripts and the chronological protocol are in
``pipeline-simulation-benchmarking/chenxi_experiments/26-capped-distance-outliers``.
Data, diagnostic traces, rejected candidate snapshots and measured outputs are
outside the tool, under
``/home/cnie/scratch_folder/results_treebench_development/medicc_distance_outliers_20261015``.
All investigation jobs used ``preempt-schwarzr``.
