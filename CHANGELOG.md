# Version 1

## Unreleased

* If FastME cannot be used (executable not found, fewer than 4 samples, or FastME fails), MEDICC2 logs a warning
  and uses neighbour joining instead of stopping
* `_branch_lengths.tsv` leaves out the normal sample given with `--normal-name` (before, only a sample named
  `diploid` was left out)
* Tree topologies are inferred with FastME (balanced minimum evolution with FastME's NNI and SPR, fixed seed) by
  default. New options `--tree-method {fastme,nj}` (`nj`: neighbour joining as before, identical results) and
  `--fastme PATH` (default `fastme` on `PATH`); in python `tree_method=` and `fastme=` in `medicc.main` and
  `run_bootstrap`, `method=` and `fastme=` in `infer_tree_topology`. Bootstrap runs use the selected method and
  `--nni-mode` starts from its tree. FastME (bioconda) is added to `doc/medicc2.yml`. An unknown method fails before the pairwise
  distances are calculated. FastME gets temporary sample names, so any sample
  name works (also long names and names with whitespace). FastME is not needed with `--tree` (unless bootstrapping), with `--tree-method nj` or for up to 3
  samples (one possible topology). Internal nodes keep their names (`internal1`, `internal2`, ...)
* Optional NNI topology refinement (`--nni-mode`, `--nni-max-iter`; `medicc.main(..., nni=True, nni_max_iter=...)`):
  steepest descent over nearest-neighbour interchanges, each candidate scored by the total tree length after
  ancestor reconstruction (`medicc/nni.py`). Writes `_nni_trace.tsv`. Works with the native model and `--fst`;
  not available with `--topology-only` and not used by bootstrap runs. Without `--nni-mode` the results are unchanged
* The native model is now the default and only built-in model (`--native` and `read_fst(native=...)` were removed).
  `medicc.io.read_fst()` returns the native model (`no_wgd`, `total_copy_numbers`); `--total-copy-numbers` and
  `--no-wgd` select its variants. `--fst PATH` / `read_fst(user_fst=PATH)` still loads a custom FST with fstlib
  (the 1.x FSTs stay in `medicc/objects`, e.g. `wgd_asymm.fst`); `--fst` with `--events` is rejected.
  `--prune-weight` > 0 is rejected unless `--fst` is given (work in progress); `--wgd-x2` has no effect.
  Removed `--input-chr-separator` and `--fst-chr-separator` (`X` is fixed) and the `chr_separator` / `wgd_x2`
  arguments of `medicc.main`. Removed the FST-based event calling and WGD detection
  (`calculate_cn_events_per_branch`, `io.load_main_fsts`, `core.detect_wgd`); `calculate_all_cn_events` now
  computes the native minimum-cost event histories
* Bugfixes
    * `_events_overlap.tsv` lists every WGD event (previously only the last WGD was kept)
    * `--bootstrap-nr` / `run_bootstrap(N_bootstrap=...)` below 5 fail early with a clear error
      (previously a division by zero during the runs)
* Native OpenFST model (opt-in with `--native` or `medicc.io.read_fst(native=True)`): losses incl. LOH and gains
  around up to two capped whole-genome doublings per branch (CN cap 8, no pruning), computed by a C++ engine
  (`medicc/native`, built by `setup.py`). Supports allele-specific and total copy numbers, `--no-wgd`, `--events`
  with BED overlaps, `--fst`, `--tree`, `--topology-only`, `-j` (threads), bootstrapping and plots.
  `--prune-weight` > 0 and chromosome separators other than `X` are rejected with `--native`; `--wgd-x2` has no
  effect. Without `--native` the results are unchanged
* Tech stack update: Python 3.12, numpy 2, pandas 3, Cython 3 and OpenFST 1.8.4
    * Packaging moved to `pyproject.toml`; `setup.py` only builds the Cython extensions.
      The upper pins `numpy<2` and `pandas<2.1` are removed. New optional extras
      `parallelization` (joblib) and `events` (pyranges)
    * `doc/medicc2.yml` (used by the GitHub test workflow) pins Python 3.12, numpy 2.5 and pandas 3.0.
      The workflow installs pyranges through the `events` extra
    * Replaced removed numpy/pandas APIs (`np.row_stack`, `np.in1d`, `read_csv(delim_whitespace=...)`,
      `DataFrame.iteritems`, positional `Series[int]` access, `int()` on a Series)
    * Explicit `observed=False` for groupings by categorical chromosomes (keeps the previous results)
    * `DataFrame.stack` no longer drops missing values in pandas 3; the input readers drop them
      explicitly as before (FASTA input with a different number of chromosomes per file, samples
      with missing segments)
    * Event tables and region overlaps are built without growing empty DataFrames, so they work with
      pandas copy-on-write and strict dtype setting. Event start, end and `cn_child` are now always
      written as integers (previously some rows were written as floats, e.g. `0.0`)
    * Removed the global `np.long` workaround for old pyranges versions
    * Plotting uses the matplotlib layout engine API instead of `set_constrained_layout_pads`

## Version 1.1.2

* Exchange `@cache` for `@lru_cache(maxsize=None)`

## Version 1.1.1

* Hotfix: Added extra compilation arguments to fix installation bugs

## Version 1.1.0

* Faster pairwise MED distance calculation 

## Version 1.0.4

* added legacy_loh_asymm_fst
* Bugfixes
    * plot_cn_heatmap in case max cn is 1
    * Bugfix for inputs with a single chromosome in _read_tsv_as_dataframe
* Misc
    * Additional checks for import_tree
    * Some changes to prevent pandas warnings

## Version 1.0.3

* Added CHANGELOG
* Notebooks
    * Bugfix in phasing_demo
    * Added "Detecting multiple WGDs" to WGD_detection notebook
* Added --no-plot-tree flag
* Misc. warnings and f-string formats

