Faster capped-WGD ancestral joins
================================

The change
----------

The expensive step is the bottom-up join of two retained child languages.
Previously its initial upper bound came from the exact pair distance between
one chosen word from each language. Those words are feasible, but other words
in the languages may admit a substantially cheaper shared ancestor. A loose
initial bound can allocate millions of states that become irrelevant once the
search finds that cheaper history.

Before searching each fixed-WGD branch pair, we now walk one feasible path
through the actual primitive transitions and both child languages. For each
ancestor symbol, the two sides are independent; we choose their cheapest
available edges according to edge cost plus the existing suffix lower bound.
We follow the best such common symbol and repeat to the end. This is a greedy
witness, not an approximate replacement for the exact search.

In Python-style pseudocode::

    upper = cost_of_existing_feasible_history
    for left_wgd, right_wgd in branch_order:
        offset = left.offset + right.offset + left_wgd + right_wgd
        witness = greedy_complete_path(left, right, upper - offset)
        upper = min(upper, offset + witness)
        result = exact_all_tie_search(left, right, upper)
        upper = min(upper, result.cost)
    return union_of_results_with_cost(upper)

The helper returns the supplied bound unchanged if it cannot complete a path.
For a completed path it includes all child arc/final weights and primitive
weights. Child offsets and the global WGD charges are added exactly once by the
caller. X is an ordinary explicit separator in this walk; the helper does not
factor a message into independent chromosomes or lose correlations between them.

Why ties survive
----------------

A completed witness is an actual permitted history, so its cost cannot be less
than the optimum. Replacing an upper bound with its smaller feasible cost cannot
exclude an optimal path. All ancestral cutoffs still use ``cost > upper``;
paths tying the bound remain eligible. The exact forward/backward search and
OpenFST determinization/minimization are unchanged. All WGD branch pairs tying
the best score are still united. Concurrent branch workers share the feasible
bound using the existing atomic incumbent.

This preserves the existing local prune_weight=0 rule. It does not change the
model into a globally unpruned optimization of all ancestors on the entire tree.
The top-down conditional reconstruction algorithm is unchanged.

Reading the code
----------------

* ``medicc/full_fst/model.py``: constructs the initial feasible join bound.
* ``medicc/cext/full_fst/common/src/feasible_path.hpp``: shared greedy walk.
  The C++ template accepts either engine's Side object, much like a Python
  function accepting either object with start/prepare/edges/final methods.
* ``medicc/cext/full_fst/ancestors/src/nfa_engine.hpp``: calls the helper before
  each exact sibling-join search; logs the initial and tightened bounds.
* ``medicc/cext/full_fst/direct_pair/src/pairwise.cpp``: reuses the same helper.
* ``tests/check_ancestor_bounds.py``: independently checks complete weighted
  languages, including ties, offsets, separators and parallel branch searches.

No new dependency or replacement graph engine is introduced. OpenFST continues
to handle the retained-language operations; the existing Boost.Graph-backed
search and custom lazy transition generation remain in place.

Validation and reproduction
---------------------------

The scratch candidate passed the literal FST oracle (including 240 ancestral
messages, incremental NNI and serial/parallel NNI), high-CN regressions, Python
API/ownership checks and 14,745 fixed-WGD pair witness/cutoff checks. An additional
192 weighted all-tie joins matched the literal oracle with one and two workers.
These weighted cases put costs on transitions and final states and use nonzero
child offsets.

Scripts live outside this tool repository in
``pipeline-simulation-benchmarking/chenxi_experiments/22-capped-ancestral-optimization``.
Inputs are referenced from experiment 21's frozen selection; no dataset is
copied into this tool checkout. Scratch snapshots, logs and complete comparisons
are under
``/home/cnie/scratch_folder/results_treebench_development/medicc_capped_ancestral_20261015``.
Every task job uses ``preempt-schwarzr``. Baseline source is preserved at commit
``672411e`` for the pair-optimized implementation and ``0c66cd9`` for the original
capped implementation.

Fixed-topology runs use eight workers and retain the original reference_pair
witness backend to isolate the ancestral-engine change. Their profiles and every
branch cost are compared against the original runs. The FastME workflow uses
the optimized direct_pair backend and includes all distances, topology inference,
ancestral reconstruction and branch scoring. The nine selected trees were chosen
as difficult cases; they are not an unbiased sample of simulation performance.

Benchmark results
-----------------

All nine fixed-topology runs retained identical full profiles and branch costs.
Ancestral reconstruction was 1.04--2.50x faster,
with a median speedup of 1.83x. Eight workers were used throughout.

.. list-table:: Fixed-topology reconstruction, seconds
   :header-rows: 1

   * - Dataset / iteration
     - Before
     - After
     - Speedup
   * - lgg_high_wgd / 9
     - 304.0
     - 154.4
     - 1.97x
   * - meso_low_wgd / 12
     - 188.6
     - 105.0
     - 1.80x
   * - meso_low_wgd / 17
     - 502.8
     - 200.9
     - 2.50x
   * - meso_low_wgd / 7
     - 1124.5
     - 559.6
     - 2.01x
   * - meso_low_wgd / 9
     - 196.1
     - 107.1
     - 1.83x
   * - pancancer / 18
     - 157.9
     - 95.1
     - 1.66x
   * - pancancer / 2
     - 227.8
     - 218.1
     - 1.04x
   * - pancancer / 9
     - 150.7
     - 108.7
     - 1.39x
   * - ucec_mid_wgd / 8
     - 2441.2
     - 1162.4
     - 2.10x

The main table uses the first completed trial per dataset; the controlled
repeat below is reported separately. Peak RSS fell on every tree. In the
largest case (ucec iteration 8), it fell from about 29.1 to 9.5 GiB.

The completed lgg_high_wgd iteration 9 full workflow (250 leaves) changed from
4,492.0 to 246.4 seconds (18.23x). Distance computation changed from 4,162.7 to
82.0 seconds; ancestral reconstruction changed from 312.8 to 147.2 seconds
(2.12x). All 31,375 distances, the topology file, all reconstructed profiles and
branch lengths were identical. RF and ancestral accuracy therefore did not
change. These wall times use separate cluster allocations; native state counts
and a separate same-allocation trial are recorded in the external report.

The nine-tree fixed-topology comparison is recorded in
``comparison.json`` / ``report.txt`` in the scratch results directory. The
same-allocation comparison is recorded in ``paired_comparison.json``. On meso
iteration 12, that sequential before/after trial changed ancestral reconstruction
from 209.3 to 92.5 seconds (2.26x) with identical profiles and branch costs. Running
trials are omitted from those reports until their output comparisons succeed.

The warm start reduces states allocated above the eventual optimum; it does not
remove the cost of exploring histories that still tie the optimum. Consequently
large ambiguous joins remain the principal limitation, and speedups vary by tree.
