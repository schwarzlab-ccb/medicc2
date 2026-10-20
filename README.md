> **Codex-generated experimental implementation.** This reference branch, including
> the capped event caller and BED-annotation integration, contains code generated
> by Codex. The checks documented below do not replace human code review.

Experimental distance optimization: this branch restores the original
short-chromosome-first policy. See [the review note](doc/pair_distance_outliers.rst)
for its provenance, exactness argument and matrix validation.

# Optimized capped-WGD pilot branch

Current branch: `full-fst-lloh-wgd2-capped-x2`. The optimization work was based on
`full-fst-lloh-wgd2-capped-x2` at `0c66cd9` (the capped model originated from
`full-fst-lloh-wgd2` at `f050793`). Mandatory WGD maps `c -> min(2*c, 8)` for CN 0–8;
X remains X. WGD costs one per genome, with at most two per directed edge.
Other events still operate on the bounded 0–8 alphabet; this does not model
hidden CN above 8.

The current implementation uses an **OpenFST-reduced AD factor**, a
**shared-ancestor pattern bound**, and faster exact suffix-table preparation.
The existing Boost.Graph A* and backward language reconstruction remain in use.
All optimal ancestor ties are retained. Downward reconstruction keeps the
original primitive order so its historical choice among tied profiles is stable.

Across nine selected difficult 250-leaf trees, median reconstruction time is
**26.4 seconds**, a median **6.06× speedup** over the earlier warm-start capped
implementation. All profiles and branch costs match. The hardest case changed
from **1,162 to 55.3 seconds**, close to its old uncapped 47.7 seconds; peak RAM
fell from about 9.5 to 2.4 GiB. A same-allocation comparison against the immediate
parent changed 107.6 to 25.3 seconds. These are selected-case measurements;
the old uncapped implementation has not received the new optimizations.

Read [the search-explosion review](doc/search_explosion_optimization.rst) for
Python-style pseudocode, the exactness proof, complete timings and limitations.
Earlier changes are documented in the [pair warm-start note](doc/pair_search_optimization.rst),
[ancestral warm-start note](doc/ancestral_search_optimization.rst), and
[forward-relaxation note](doc/search_relaxation_optimization.rst).

Model assets use the same capped mapping throughout. Regenerate chromosome
bounds with `python scripts/generate_chromosome_bound.py`; generate the compact
factor with `python scripts/generate_capped_ad.py --help`, then build as below.
The independent OpenFST oracle uses `wgd_cap=True`; legacy factory calls keep
their previous default. Pilot scripts, datasets and results remain outside this
tool checkout. The inherited PDF explains the parent model and predates these
optimizations; use the review notes for the current capped implementation.

# MEDICC full L_LOH / two-WGD reference branch

This is a standalone MEDICC source tree for reviewing the full-FST implementation.
The normal entry point is `medicc2 INPUT.tsv OUTPUT_DIR`. It accepts ordinary
MEDICC sample names and input tables, adds a diploid normal when needed, and uses
FastME followed by full-model ancestral reconstruction. `--tree` accepts a supplied
topology; `--nni-mode` enables MEDICC's incremental NNI search.

The original `full-fst-lloh-wgd2` reference branch of `medicc2` was based on
`imbalance_FST_tree_search_dev` at `c9dddb8c36d9`. The original working
checkout and the `medicc2-pre-pub` checkout remain unchanged. The copied core and ancestor
modules now dispatch directly to `medicc.full_fst.ancestors`; runtime monkeypatch
installation is no longer needed. The full model is the default of this branch's
CLI. This is a reference implementation, not yet a complete replacement for every
legacy CLI option. Custom FSTs, event-history output, bootstrapping, length-aware
selection and total-CN mode are not exposed by this CLI.

The branch imports the reviewed standalone source tree. Legacy example/output
datasets and the former FSTlib submodule arrangement remain in the parent Git
history; this branch vendors the FSTlib source and its license.

## Build and run

The existing environment works on this machine:

```bash
MEDICC_PYTHON=/projects/ag-schwarzr/chenxi_conda/medicc2_pre_pub_env/bin/python
$MEDICC_PYTHON -B scripts/build.py --legacy
$MEDICC_PYTHON -B medicc2 /path/to/input.tsv /path/to/output \
  --normal-name internal_0 -j 8 \
  --fastme /projects/ag-schwarzr/chenxi_conda/medicc3_dev_env/bin/fastme \
  --no-plot-tree
```

`--legacy` rebuilds the copied FSTlib and MEDICC Cython extensions; the three new
native engines are built in all ordinary builds. Builds need a C++17 compiler,
Python, Cython, OpenFST headers/libraries, and Boost.Graph headers (tested with
OpenFST 1.8.4 and Boost 1.75). Boost.Graph is header-only for these operations;
it adds no runtime library. The explicit conda environment below does not include
the system Boost headers used on this machine; provide Boost headers separately
when building on another machine. `environment-linux-64.explicit.txt`
records the original Linux conda environment plus FastME 2.1.6.3:

```bash
conda create -p /path/to/new/env --file environment-linux-64.explicit.txt
conda activate /path/to/new/env
python -B scripts/build.py --legacy
```

Optional installation into that **separate environment**:
`python -m pip install --no-build-isolation .`. It installs a `medicc2` command and
builds the native libraries. Running directly from this source tree requires no
installation and does not replace another MEDICC package.

Typical options:

```bash
python -B medicc2 INPUT.tsv OUT --normal-name NORMAL --tree TOPOLOGY.nwk -j 8 --no-plot-tree
python -B medicc2 INPUT.tsv OUT --normal-name NORMAL --nni-mode -j 8 --fastme FASTME --no-plot-tree
```

Output uses MEDICC names: `<prefix>_pairwise_distances.tsv`,
`<prefix>_final_tree.new`, `<prefix>_final_cn_profiles.tsv`,
`<prefix>_branch_lengths.tsv`, and `<prefix>_summary.json`.
No distance matrix is computed for a supplied tree. Ancestor profiles and NNI
steps are also saved. All inputs and output directories are caller supplied.

Python API:

```python
import medicc
from medicc.full_fst import reconstruct, PairScorer

frame = medicc.io.read_and_parse_input_data("INPUT.tsv", normal_name="NORMAL")
result = reconstruct(frame, normal_name="NORMAL", n_cores=8, fastme="fastme")
# result.tree, result.copy_numbers, result.pairwise_distances
with PairScorer() as scorer:
    distance = scorer.score("1X1", "8X8")
    directed = scorer.score("1X1", "8X8", directed=True)
```

## Review entry points

Read this README, then [the Python-first review guide](doc/full_fst_review.rst),
then [the algorithm slides (PDF, 22 slides)](doc/full_fst_algorithms.pdf).
The guide covers the module map, ownership rules and a C++ translation key;
the slides explain the model, searches and division of work between Python,
custom C++, OpenFST and Boost.Graph. Then start the code review at
[`reconstruct()` in medicc/full_fst/core.py](medicc/full_fst/core.py).
For the command-line entry, follow `medicc2` into `medicc/full_fst/cli.py`.


- `medicc2` and `medicc/full_fst/cli.py`: general tool entry point.
- `medicc/full_fst/core.py`: input encoding, distances, FastME and reconstruction.
- `medicc/full_fst/ancestors.py`: the former NNI hotpatch, now an explicit backend.
- `medicc/core.py` and `medicc/ancestors.py`: small dispatch additions to the frozen MEDICC code.
- `medicc/full_fst/model.py`, `message.py`, `scheduler.py`: native message ownership and tree passes.
- `medicc/cext/full_fst/direct_pair/src/pairwise.cpp`: exact factorized pair scorer.
- `medicc/cext/full_fst/ancestors/src/nfa_engine.hpp` and `ancestors.cpp`: all-minimum joins and top-down traceback.
- `medicc/full_fst/_bindings.py`: centralized ctypes signatures and library loading.
- `source_provenance.json`: historical source paths and hashes before subsequent adaptations; current hashes are in `code_manifest.json`.

## Native algorithm ownership

The production C++ kernels use OpenFST for determinization, weighted union,
and epsilon removal. The measured acyclic integer minimizer and weight
normalization remain specialized: the fully native combination exceeded the
5% runtime budget in paired tests. Boost.Graph's
`astar_search_no_init` handles ancestral-message search and conditional
ancestor traceback. Ancestral state indexing retains the specialized index after measuring
standard hash-container alternatives. The
fixed-profile pair scorer retains its specialized search and state index: the
generic graph adaptation exceeded the runtime budget for distance-matrix
construction followed by FastME.

MEDICC supplies the model, admissible bounds, bounded lazy composition, WGD and
chromosome bookkeeping, and extraction of the complete tight-path ancestor
language. The composition engine remains specialized because the tested native
OpenFST composition alternatives exceeded the performance budget.

`medicc/cext/full_fst/common/src/openfst_language.cpp` adapts the integer
acceptors to OpenFST and is compiled separately from the search kernels.
`boost_search.hpp` adapts the implicit bounded graph and preserves the priority
ordering used by the tool. The Python bindings remain `ctypes` over a C ABI;
legacy FSTlib still uses its Cython/OpenFST wrappers. `reference_pair` remains an
independent custom implementation for comparisons.

The upstream algorithm implementations are [OpenFST](https://www.openfst.org/)
and [Boost.Graph](https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/astar_search.html).
The shared headers and translation unit here adapt their APIs to MEDICC.

Library evaluation and reproduction scripts live separately in
`pipeline-simulation-benchmarking/chenxi_experiments/20-reproduce-medicc-figures/library_reuse`.
No experiment inputs or results are included in this tool package.

## Model and retained experiment semantics

Let `A = L_LOH * G`. A directed edge allows `A`, `A * D * A`, or
`A * D * A * D * A`, where `D` is a mandatory doubling with cost one.
The limit is **two WGDs per directed edge**, not two per whole tree. CN states
are 0–8; input CNs above `--maxcn` are clipped by MEDICC's existing reader.
Chromosome and allele boundaries are `X` separators. The pair kernel minimizes
over a shared ancestor, with global WGD counts fixed across all chromosomes on
each side before taking their minimum. It does not choose independent WGD counts
per chromosome.

Pair and branch scoring may collapse identical adjacent column pairs. Ancestor
messages retain the complete input bin grid and all tied local minima at
`prune_weight=0`. This preserves the experiment's local reconstruction semantics;
it does not claim an unpruned global optimum over all ancestors of the tree.
No arbitrary native search-state or input-length cutoffs are introduced.

FastME uses BME, NNI and SPR, seed 104000. MEDICC NNI is opt-in, defaults to one
expanded frontier tree, and stops on a plateau; its evaluation limit is unlimited
unless explicitly supplied. The factorized scorer is the default. The retained
`reference_pair` implementation is used by experiment 20 to reproduce the original
truth-tree workflow exactly.

## Validation and figure reproduction

Tool-only checks (no input datasets are included):

```bash
python -B tests/check_pairs.py
python -B tests/check_oracle.py
python -B tests/check_ancestor_bounds.py
python -B tests/check_api.py
python -B tests/check_bindings.py
```

The 900-dataset definitions, run-array script, external-reference checks, and
collectors belong in the separate
`pipeline-simulation-benchmarking/chenxi_experiments/20-reproduce-medicc-figures`
folder. This repository contains no experiment input datasets or result files.
The small `ancestral_parts.txt` files are fixed model transition tables needed by
the native implementation, not experiment data.

Validation on 14 October 2026: 52,488 exhaustive directed/kernel comparisons,
240 ancestral-message checks, 64 incremental-NNI neighbors, serial/parallel NNI,
general sample-name and normal-creation checks, nine historical workflow runs,
and three actual CLI runs passed. The wheel built and passed an isolated-install
native smoke check. The complete 900-task arrays were not rerun.

External validation reports are under
`/home/cnie/scratch_folder/results_treebench_development/medicc_full_fst_reference_review_20261014/`:
`check_1791955321920609818_report.json` and `cli_1791955335515184876_report.json`.

The library-backed ancestral implementation was checked against the preserved
pre-replacement reference on nine workflows, with five paired trials each.
All 90 executions matched full profiles, trees/branch lengths, scores and NNI
traces. Median reconstruction runtime changes were +2.05% to +2.81% on this host
(single worker; input/import time excluded). Synthetic weighted-language,
independent inference, API and concurrency checks also passed. This is measured
workload coverage, not a guarantee for every dataset or machine.

The separate library evaluation report, raw timing/memory records and rejected
alternatives are under
`/home/cnie/scratch_folder/results_treebench_development/medicc_openfst_optimization_20261014/`:
`library_reuse_assessment.txt` and
`library_reuse_acceptance_reuse_candidate_index.json`.
The source-only replay patch belongs to experiment 20's `library_reuse` scripts.

The readability refactor is documented in `doc/full_fst_review.rst`. Its source
patch and paired-workflow validation scripts are in experiment 20's
`review_refactor` folder. Current before/after checks and wheel validation are
recorded separately under
`/home/cnie/scratch_folder/results_treebench_development/medicc_readability_review_20261014/`.
The earlier benchmark reports above describe their preserved implementation
snapshots; their source paths have not been retroactively rewritten.

## Optional capped-model event calling

Add `--events` to the reference CLI (or `events=True` to `reconstruct`) to call
one minimum-cost event history for each reconstructed branch. This uses the
same capped full L_LOH model, including up to two genome-wide WGDs; it does not
load the legacy event-calling FST. Every history is replayed and its event count
must equal the directed branch score. Tied histories prefer fewer WGDs.

Existing output filenames stay unchanged. Events are written to MEDICC's usual
`<prefix>_copynumber_events_df.tsv`, with the usual leading columns and additional
order/stage fields. See `doc/full_fst_review.rst` for their meaning. Calling events
adds work after reconstruction and does not change the tree or ancestor profiles.
Run `python -B tests/check_events.py` for exhaustive short-genome and workflow checks.

Event-caller validation (20 October 2026): 13,122 exhaustive two-bin histories
matched the independent literal-FST oracle; 52,488 existing pair-score checks
and API/NNI/CLI regressions passed. Saved 20- and 100-leaf lambda-5 trees produced
173 and 972 replayable unit events, exactly matching their branch scores. The
20-leaf full CLI also reproduced the saved profiles and branch lengths. Logs and
benchmark event tables are outside the tool in
`/home/cnie/scratch_folder/results_treebench_development/medicc_capped_events_20261020/`.

### Annotating events with chromosome arms and genes

With `--events`, add `--chromosomes-bed FILE` and/or `--regions-bed FILE`.
Use `default` for the bundled hg38 chromosome arms or Davoli gene list; `none`
(the default) disables that annotation. For example, from this checkout:

```bash
python medicc2 input.tsv output --events \
  --chromosomes-bed medicc/objects/hg38_chromosome_arms.bed \
  --regions-bed medicc/objects/Davoli_2013_TSG_OG_genes.bed
```

The three BED files in `medicc/objects` are unchanged copies from MEDICC:
`hg19_chromosome_arms.bed`, `hg38_chromosome_arms.bed`, and
`Davoli_2013_TSG_OG_genes.bed`. Choose files matching the input's genome build;
the command does not convert coordinates. BED input requires four tab-separated
columns: chromosome, start, end, and region name, without a header (`#` comments
are allowed).

This writes the usual `<prefix>_events_overlap.tsv` after calling capped-model
events. It labels regions with strictly more than 90% coverage by gains or losses,
separately for each branch and allele. It does not change reconstruction or the
raw event table. WGDs remain separate rows, including two on the same branch.

Annotations use MEDICC's existing PyRanges overlap code. Install the optional
package extra (`pip install '.[events]'`) or, for an already built checkout,
`python -m pip install 'pyranges>=0.1,<0.2'`. Event calling without BED annotations
does not require PyRanges. Validate with `python -B tests/check_event_annotations.py`.


### FastME taxon-name compatibility

Names at or above64UTF-8bytes are passed to FastME as unique short IDs. The original names are restored before MEDICC rooting and ancestral reconstruction. `fastme_taxa.tsv` records the bijection; `fastme_short_ids.nwk` retains FastME's short-ID output; `fastme_raw.nwk` contains the same tree with original names. Shorter supported names retain the existing interface. The64byte boundary is mapped as well because the installed FastME build produced invalid names there. Names are never truncated.

Regression command: `python -s tests/check_fastme_labels.py` (optionally set `FASTME` to the binary path). Tests cover63/64/65byte names, colliding prefixes, original names resembling short IDs, restored branch distances, mismatched matrix axes, duplicate inputs and unexpected output taxa.
