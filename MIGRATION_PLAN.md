# Migration plan: capped full-FST (L_LOH, ≤2 WGD/edge, CN cap 8) onto master

| | Path | Branch | Conda env |
|---|---|---|---|
| **Left monitor (REF)** | `~/Projects/project-medicc/medicc-full-fst-lloh-wgd2-capped-x2` | `full-fst-lloh-wgd2-capped-x2` (codex) | `medicc-codex-ref` |
| **Right monitor (NEW)** | `~/Projects/project-medicc/feature-full-fst-capped-wgd2` | `feature-full-fst-capped-wgd2` (from `master` @ `04e9bfe`) | `medicc-migration` |

This file is untracked. Commit it or leave it out, as you prefer.

---

## Part A: what the codex branch adds

`master` is MEDICC2 1.1.2. It uses OpenFST transducer files (`objects/*.fst`), NJ trees, and FSA intersection for
ancestors. The codex branch (74 commits ahead) contains **two layers of work**:

1. **Your own pre-codex work** (up to `c9dddb8`, from `imbalance_FST_tree_search_dev`): SPR/MCMC, NNI tree search,
   length-encoding FSTs, constrained Sankoff, and the **separate up-pass/down-pass FSTs** (`upper_pass_fst` /
   `lower_pass_fst`).
2. **Codex work** (`f050793` → `9ddbcf1`): a new, separate engine called **`medicc.full_fst`**. This is the new feature:

| Capability | Where (REF) |
|---|---|
| **Model**: `A = L_LOH·G`; per directed edge `A`, `A·D·A` or `A·D·A·D·A`, where `D` = mandatory WGD with `c → min(2c, 8)` (capped), cost 1 per genome. CN alphabet 0–8, `X` = chromosome/allele separator. | `common/src/wgd_model.hpp`, `*/assets/ancestral_parts.txt` |
| **Exact pair distance** `K(x,y)=min_z T(z,x)+T(z,y)` and **directed** score `T(p,c)`, factorised per chromosome with *global* WGD counts | `cext/full_fst/direct_pair/` + `full_fst/pair_scorer.py` |
| **Bottom-up join**: returns *all* locally optimal ancestors as a weighted DFA ("message", prune_weight = 0), using A* (Boost.Graph), admissible bounds, a sibling pattern bound, and segmented fallback | `cext/full_fst/ancestors/` + `full_fst/model.py`, `message.py` |
| **Top-down resolve**: picks one child profile given a fixed parent | `ancestors/src/ancestors.cpp` |
| **Tree pass scheduler** (thread-parallel, ready-queue) | `full_fst/scheduler.py` |
| **Workflow**: encode → distances → tree → ancestors → branch lengths; FastME option | `full_fst/core.py` |
| **Incremental NNI** on the native messages (reuses unchanged subtree messages) | `full_fst/ancestors.py` + `medicc/core.py::nni_mode` + `medicc/nni.py` |
| **Exact event calling** (one minimum-cost history, replay-verified) | `direct_pair/src/event_trace.hpp` + `full_fst/events.py` |
| **CLI** `medicc2` → full model only | `full_fst/cli.py` |
| Asset generators + independent literal-FST test oracle | `scripts/generate_*.py`, `tests/oracle.py`, `tests/check_*.py` |

The codex branch also **deletes or rewrites a lot**: examples, notebooks, CI, the `.fst` objects, `test_medicc.py`, and
the git submodule (it vendors `fstlib`). It also rewrites many log messages and changes CLI defaults. **We do not
migrate any of that.**

### What we migrate, and what we don't

**Migrate:** `medicc/full_fst/*` (most of it), the native `ancestors` and `direct_pair` engines, assets, two generator
scripts, the oracle, topology-only NNI helpers, and a few small additions to legacy files.

**Skip:**
* `reference_pair/`: an independent copy for experiment 20, not used in production. Optional appendix.
* `*/src/pair_engine.hpp`: dead code; nothing includes it.
* vendored `fstlib/`: byte-identical to master's submodule `8216876`. Keep the submodule.
* `medicc/cext/acyclic_prune.*`, `medicc/fst_utils.py`: these speed up the *legacy* OpenFST path only.
* `full_fst/provenance.py`, `integration.py`, `pipeline.py`: hash manifests and old-import shims.
* SPR, length encoding, constrained Sankoff, `la_nni`, and the up/down FST split.
* all cosmetic edits to `io.py`, `plot.py`, `tools.py`, `sim.py`, `stats.py`, `bootstrap.py`.

---

## Part B: rules for this migration

### B1. Keep master's single FST for both passes (drop the upper/lower split)
master's `reconstruct_ancestors(tree, samples_dict, fst, normal_name, ...)` already uses **one** FST for the up and
down passes. The `upper_pass_fst` / `lower_pass_fst` split came from your dev branch (`5f5b906`) and survives in REF.
In REF the up/down split survives only as the sentinel `FULL_TWO` passed **twice** (`upper_pass_fst=FULL_TWO,
lower_pass_fst=FULL_TWO`). That keeps the old legacy signature alive. The native engine already uses one model:
`join` (up), `resolve` (down) and branch scoring all load the same `ancestral_parts.txt`. In NEW:
* **No `FULL_TWO`, no `upper_pass_fst` / `lower_pass_fst` / `fst` arguments, no `la_nni`** anywhere in `full_fst`.
* **Leave `medicc/core.py` and `medicc/ancestors.py` exactly as on master.** No dispatch hooks into legacy code.
  The full-FST stack calls its own functions directly.
* `nni_mode` and the neighbour-evaluation functions move into a new file `medicc/full_fst/nni_search.py`.
  `medicc/nni.py` keeps only pure tree-topology code.

### B2. Master's behaviour wins
Wherever a file exists on master, keep master's behaviour (log text, logger names, defaults, output files). For the new
CLI, keep master's conventions wherever an option or output also exists in master's `medicc2`:

| Behaviour | master | codex | **NEW** |
|---|---|---|---|
| `medicc2` command | legacy OpenFST | full model only | **full model is the new default `medicc2`**; master's script moves, unchanged, to `medicc2-legacy` (Step 44) |
| Logging | `logging_conf.yaml` format; `-v` / `-vv` / `--silent` set the `medicc` logger | `logging.basicConfig` (different format, root handler), `-v` only | **master** |
| `--plot` default | `auto` | `none` | **`auto`** |
| `-j/--n-cores` default | `None` | `1` | **`None`** (mapped to 1 internally) |
| Existing output prefix | overwrite | refuse | **overwrite** |
| Summary file | `_summary.tsv` (`summarize_patient`) | `_summary.json` | **`_summary.tsv`**; extra json optional, under another name |
| Pairwise file with `--tree` | written (zeros) | not written | **written (zeros)** |
| Tree inference | NJ (`medicc.nj`) | FastME BME/NNI/SPR, seed 104000 | **NJ default**; FastME opt-in via `--tree-method fastme` (decision D1) |
| Internal node names | `internalN` (from `medicc/nj.py`) | `inferred_N` | **`internalN`** (master's `--tree` check relies on the substring `internal`) |
| Output tree shape | unnamed root → (MRCA, `diploid` leaf), MRCA first | root *named* `diploid` with one child | **master shape** (wrap before writing) |
| Sample order | input order | normal first, then natural sort | **input order** |
| Final message | `logger.info("MEDICC2 has finished.")` | `print("MEDICC full-FST finished…")` | **master** |
| Launcher | nothing | sets `OPENBLAS/OMP/MKL_NUM_THREADS=1`, `dont_write_bytecode` | **nothing** |
| `tools.py` logger | `'medicc.io'` | `__name__` | **master** (only add `sum_of_branch_length`) |

### B3. Typing vs copying
Type all hand-written code. **Copy, don't type** generated tables (`generated/*.hpp`, `ancestral_parts.txt`, ≈60k
lines). They are model data, and Steps 4–5 regenerate them to prove the copies are correct. When a file is
byte-identical in `ancestors/` and `direct_pair/`, type it once and `cp` it (the step says so and checks with `cmp`).

### Decisions I made for you (change them if you disagree)
* **D1.** NJ stays the default tree method (master behaviour); FastME is opt-in.
* **D2.** NNI CLI defaults: master has no NNI, so I kept codex's (`--nni-max-topology` unlimited, `--nni-max-frontier 1`,
  stop on plateau). Your pre-codex branch used `--nni-max-iter 1000`.
* **D3. (decided: fix it)** Step 38 ports codex's fix to a real master bug in `event_reconstruction.overlap_events`
  (codex found it; reproduced on master's code in `medicc-migration`). Every WGD is written to the same row label
  `('all','0','0')`, so with WGDs on two branches `_events_overlap.tsv` keeps only the last one. Only that BED-overlap
  summary is affected; `_copynumber_events_df.tsv` is fine. The function is shared, so `medicc2-legacy` gets the fix
  too: its output changes only when more than one WGD exists, i.e. exactly where it was wrong. Worth reporting upstream.

* **D4. (decided: accept `0`)** Phase 0 / Step 0.5: old master (pandas 2) writes `cn_child` of *some* WGD rows as
  `0.0` and all other rows as integers. It's a pandas type-inference quirk; on the new stack with the int cast it's
  always `0`. The value is the same, so we leave it: no special-case code to mimic pandas 2, and no switch to REF's
  `2.0`-everywhere. That one cell is the **only** output difference between old- and new-stack master across every
  example run (see Step 0.7).

---

## Step 0: environments (done by Claude, verified)

The two environments are **identical** (same package list, checked with `mamba list --explicit`):

* **`medicc-codex-ref`**: REF's exact `environment-linux-64.explicit.txt` (Python 3.12.13, numpy 2.5.1, pandas 3.0.3,
  biopython 1.87, matplotlib 3.11.1, Cython 3.1.2, OpenFST 1.8.4, FastME 2.1.6.3). Plus `libboost-headers` 1.92 (REF's
  file leaves Boost out, but the ancestor engine needs it), plus `pytest` 9.2 + `pytest-timeout` (conda) and
  `pyranges` 0.1.7 (pip) for master's test suite and BED annotations. REF is built here (`scripts/build.py --legacy`),
  and its checks pass: `check_pairs`, `check_resolve_tie`, `check_api`, `check_bindings`, `check_events`.
* **`medicc-migration`**: a clone of `medicc-codex-ref` plus the same three test packages. **REF is not installed or
  built into it.** (Neither env has `medicc` pip-installed; both run from their source tree.)

The old toolchain (Python 3.10, numpy 1.26, pandas 2.0.3) is still available in your existing env
`medicc2-imbalance-fst-tree-search-dev-env`. Phase 0 uses it only as the "before" side of an output comparison.

Daily routine:
```bash
# left terminal                                   # right terminal
mamba activate medicc-codex-ref                   mamba activate medicc-migration
cd ~/Projects/project-medicc/medicc-full-fst-lloh-wgd2-capped-x2
                                                  cd ~/Projects/project-medicc/feature-full-fst-capped-wgd2
```
Playing with REF: `python -B medicc2 IN.tsv OUT --fastme $CONDA_PREFIX/bin/fastme` (use master's
`examples/simple_example/simple_example.tsv`, since REF deleted its examples). The guide is
`doc/algorithm_reading_guide.rst` (or the PDF). Read pages 1–10 before Phase 2.

Migration tests go in **`tests/migration/`** in NEW as plain scripts (pytest is only used for master's own suite). Every Python test
starts with this header (written once here, referred to as *the header* below):
```python
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
```
Run tests from NEW's root: `python -B tests/migration/<file>.py`. Commit after each green step.

---

# Phase 0: move master onto the new toolchain (no full-FST code yet)
**Goal.** master, *unchanged in behaviour*, running on REF's stack: Python 3.12, numpy 2.5, pandas 3.0, Cython 3.1,
OpenFST 1.8.4, Biopython 1.87. Every later phase builds on this. Only **toolchain** fixes go in. Codex's legacy-file
rewording, logger renames and the `future_stack=True` additions (already pandas 3's default, so they do nothing) are
**not** ported.

All of this was tried on a scratch copy of master first:
* **Before** the fixes: 9 of 11 core integration tests fail, from 3 root causes.
* **After** the fixes: `test_fsts.py` 66/66, core integration 11/11, and Gundem 40/40 pass.
* Output comparison against old-stack master: see Step 0.7.

## Step 0.1: record "before" outputs on the old stack (unmodified master)
**What/why.** You need proof that the upgrade changes nothing. Generate master's outputs on the old toolchain
*before* editing anything, then compare after Step 0.6. Your dev env has the old stack but no pyranges, so make a
throwaway clone. The old (cpython-310) and new (cpython-312) fstlib builds can live side by side in `fstlib/cext/`.
```bash
mamba create -y -n medicc-oldstack --clone medicc2-imbalance-fst-tree-search-dev-env
mamba run -n medicc-oldstack python -m pip install 'pyranges>=0.1,<0.2'
mamba activate medicc-oldstack && python setup.py build_ext --inplace && mamba deactivate
```
Type the comparison script **`tests/migration/compare_stacks.sh`**:
```bash
#!/usr/bin/env bash
# usage: tests/migration/compare_stacks.sh old|new    (run from the repo root, in the matching env)
#        tests/migration/compare_stacks.sh diff
set -u
OUT=/tmp/medicc_stack_compare
EV="--events --chromosomes-bed default --regions-bed default --plot none --no-plot-tree"
runs=(
  "simple|examples/simple_example/simple_example.tsv $EV"
  "simple_j4|examples/simple_example/simple_example.tsv --n-cores 4 $EV"
  "testing|examples/testing_example/testing_example.tsv $EV"
  "testing_tcn|examples/testing_example/testing_example.tsv --total-copy-numbers --input-allele-columns cn_a $EV"
  "testing_nowgd|examples/testing_example/testing_example.tsv --no-wgd $EV"
  "testing_x2|examples/testing_example/testing_example.tsv --wgd-x2 $EV"
  "wgd|examples/WGD_example/WGD_example.tsv $EV"
  "ov|examples/OV03-04/OV03-04_descr.txt -i fasta --normal-name OV03-04_diploid $EV"
)
for p in PTX004 PTX005 PTX006 PTX007 PTX008 PTX009 PTX010 PTX011 PTX012 PTX013; do
  runs+=("g_$p|examples/gundem_et_al_2015/${p}_input_df.tsv $EV")
done
LEGACY=${LEGACY:-medicc2}          # becomes medicc2-legacy after Step 44
if [ "$1" = diff ]; then
  diff -r -x '*.png' -x '*.pdf' -x '*.log' $OUT/old $OUT/new && echo "PASS: identical outputs"; exit
fi
mkdir -p $OUT/$1
for r in "${runs[@]}"; do
  name=${r%%|*}; args=${r#*|}
  python $LEGACY $args $OUT/$1/$name > $OUT/$1/$name.log 2>&1 || echo "FAIL $name (see $OUT/$1/$name.log)"
done
echo "done: $OUT/$1"
```
**Run:** `mamba activate medicc-oldstack && bash tests/migration/compare_stacks.sh old`. Takes a few minutes; it should
print no `FAIL` lines.

## Step 0.2: see it fail on the new stack ("red")
```bash
mamba activate medicc-migration
python setup.py build_ext --inplace          # builds the cpython-312 fstlib extensions
python -m pytest medicc/test_fsts.py -q       # 66 passed
python -m pytest medicc/test_medicc.py -q -k "not gundem and not ipynb" --timeout=600
```
Expected: **9 failed, 2 passed**. The 3 root causes (visible in the subprocess tracebacks) are fixed by Steps 0.3–0.5:
`AttributeError: module 'numpy' has no attribute 'row_stack'`,
`TypeError: read_csv() got an unexpected keyword argument 'delim_whitespace'`,
`TypeError: Invalid value '[5]' for dtype 'float64'`.

## Step 0.3: `medicc/plot.py`: `np.row_stack` → `np.vstack` (removed in numpy 2.0)
**REF:** `medicc/plot.py:372`. One line. (REF's other edits in this file are message rewording plus the
NNI-only `plot_nni_trace`; skip them.)
**Test:** `python medicc2 examples/simple_example/simple_example.tsv /tmp/p03 --plot both` → creates
`simple_example_cn_profiles.pdf` and `simple_example_cn_profiles_heatmap.pdf`.

## Step 0.4: `medicc/io.py`: `delim_whitespace=True` → `sep=r"\s+"` (removed in pandas 3.0)
**REF:** `medicc/io.py:246`. **Deviation:** write a *raw* string `r"\s+"`. REF writes `'\s+'`, which Python 3.12
reports as `SyntaxWarning: invalid escape sequence '\s'` (visible in REF's test output). Skip REF's
`future_stack=True` additions and message rewording.
**Test:** `python -W error::SyntaxWarning -c "import medicc"`, then
`python medicc2 examples/OV03-04/OV03-04_descr.txt /tmp/p04 -i fasta --normal-name OV03-04_diploid --plot none`
(the FASTA reader is the only caller).

## Step 0.5: `medicc/event_reconstruction.py`: `cn_child` under pandas 3
**What/why.** `events_df` is grown with `pd.concat([events_df, pd.DataFrame(index=…)])`. In pandas 3 that makes
`cn_child` a float64 column, and assigning an object array into it raises. REF casts the values to float (line 314),
which fixes the crash but **changes master's output**: every `cn_child` is written as `2.0` instead of `2`.
**Type:** REF's `.astype(float)` on line 314, **plus one deviation line** right before
`events_df['chrom'] = tools.format_chromosomes(events_df['chrom'])`:
`events_df['cn_child'] = events_df['cn_child'].astype(int)`. This restores master's integer output.
**Test:** `python -m pytest medicc/test_medicc.py -q -k "simple_example or testing_example" --timeout=600`, then
```bash
python medicc2 examples/simple_example/simple_example.tsv /tmp/p05 --events --plot none --no-plot-tree
cut -f7 /tmp/p05/simple_example_copynumber_events_df.tsv | grep -c '\.'      # 0: integers, as on master
```

**Known difference (D4, accepted):** old master's output has one pandas 2 quirk: the `cn_child` of *some* WGD rows is written
`0.0` (it depends on the rows already in the table; e.g. `testing_example`, `WGD_example` and OV03-04 have one each,
while Gundem PTX007/PTX012 have WGDs but none). With the int cast the new stack writes `0`. Values are equal; only that
cell's text differs. We accept this; don't add code to reproduce `0.0`.

## Step 0.6: packaging and environment files
* `setup.py`: `'numpy>=2.0,<3'`, `'pandas>=3.0,<4'`, add `python_requires='>=3.12'` (keep master's layout and the
  other pins).
* `pyproject.toml`: drop `"pathlib"` from `requires`. It is a dead PyPI backport, and pip would try to install it
  into the build environment on Python 3.12.
* `doc/medicc2.yml` (used by CI): pin the tested versions: `python=3.12`, `numpy=2.5.1`, `pandas=3.0.3`,
  `cython==3.1.2`, `openfst=1.8.4`, and **`biopython=1.87`**. Without that pin the CI solve picks Biopython 1.88,
  which writes Newick branch lengths as `1` instead of `1.00000`, so CI outputs would differ from your local env.
  (Checked with `mamba env create --dry-run`.) `libboost-headers` and `fastme` are added in Step 46, when the native
  engines need them.
* `environment-linux-64.explicit.txt` (new file, like REF's): `mamba list -n medicc-migration --explicit > environment-linux-64.explicit.txt`.
  Add a comment line saying `pyranges` 0.1.7 is pip-installed on top.

**Test:** build a wheel *without installing it*:
```bash
python -m pip wheel --no-build-isolation --no-deps -w /tmp/whl . && ls /tmp/whl   # medicc2-1.1.2-cp312-…whl
rm -rf build medicc2.egg-info
```

## Step 0.7: green, plus identical outputs
```bash
python -m pytest medicc/test_fsts.py -q                                               # 66 passed
python -m pytest medicc/test_medicc.py -q -k "not gundem and not ipynb" --timeout=600 # 11 passed
python -m pytest medicc/test_medicc.py -q -k "gundem and not ipynb" --timeout=900     # 40 passed (~5 min)
bash tests/migration/compare_stacks.sh new && bash tests/migration/compare_stacks.sh diff
```
Expected `diff`: 22 of 28 runs identical. The other 6 (`testing`, `testing_tcn`, `testing_j4`, `testing_x2`, `wgd`,
`ov`) differ only in the one D4 `cn_child` cell (`0.0` old vs `0` new) in `*_copynumber_events_df.tsv`; that is
expected and accepted. Any other difference is a real regression. Then `mamba env remove -n medicc-oldstack`, and commit:
"Move master to Python 3.12 / numpy 2 / pandas 3 toolchain".

---

# Phase 1: the model in Python (no native code yet)

## Step 1: `medicc/factory.py` (edit): capped WGD in the legacy FST factory
**What/why.** The test oracle (Step 2) builds a literal one-step WGD transducer. The capped model needs
`c → min(2c, 8)` instead of "exactly ×2", so a CN of 5 can still double (to 8).
**REF:** `medicc/factory.py:116-166` (`create_1step_WGD_fst`).
**Type in NEW:** only the `wgd_cap=False` keyword, the docstring, the `if wgd_cap and not wgd_x2: raise ValueError(...)`
guard, and the two changed `cns[t] == (min(2*cns[s], max(cns.values())) if wgd_cap else 2*cns[s])` conditions.
Leave everything else in `factory.py` alone (REF's `scale` / `det_min` / whole-chrom changes are length encoding,
which we skip).
**Test:** `tests/migration/test_01_wgd_cap.py`
```python
# <the header>
import math
import fstlib
import medicc

sym = medicc.create_symbol_table(8, "X")
fsa = lambda s: fstlib.factory.from_string(s, isymbols=sym, osymbols=sym)

def score(W, a, b):
    return float(fstlib.score(W, fsa(a), fsa(b)))

capped = medicc.factory.create_1step_WGD_fst(sym, "X", wgd_cost=1, minimize=True, wgd_x2=True, wgd_cap=True)
plain = medicc.factory.create_1step_WGD_fst(sym, "X", wgd_cost=1, minimize=True, wgd_x2=True)
assert score(capped, "5", "8") == 1          # saturation: 5 -> min(10, 8)
assert math.isinf(score(plain, "5", "8"))    # master behaviour unchanged
assert score(capped, "4", "8") == 1
assert score(capped, "2X3", "4X6") == 1      # one WGD for the whole genome
assert score(capped, "3", "3") == 0          # the WGD is optional in this 1-step FST
try:
    medicc.factory.create_1step_WGD_fst(sym, "X", wgd_cap=True)
except ValueError:
    pass
else:
    raise AssertionError("wgd_cap without wgd_x2 accepted")
print("PASS step 1")
```

## Step 2: `tests/oracle.py`: the independent literal-FST oracle
**What/why.** A brute-force dynamic program over the *legacy* OpenFST primitives (`L_LOH`, `G`, `W1`). Every native
engine is later checked against it. It's 53 lines, and the best way to learn what "T(p,c)" means.
**REF:** `tests/oracle.py`. **Type exactly** into NEW `tests/oracle.py`.
**Test:** `tests/migration/test_02_oracle.py`
```python
# <the header>
import numpy as np
from oracle import get_parts, Primitive, relation, minplus

parts, sym = get_parts()
assert {"L_LOH", "G", "W1"} <= set(parts)
lab = lambda s: [sym.find(c) for c in s]
assert Primitive(parts["W1"]).score(lab("5"), lab("8")) == 1
assert Primitive(parts["L_LOH"]).score(lab("22"), lab("11")) == 1   # one interval loss
assert Primitive(parts["L_LOH"]).score(lab("2X2"), lab("1X1")) == 2 # X breaks intervals
words = list("012345678")
A = minplus(relation(Primitive(parts["L_LOH"]), words, sym), relation(Primitive(parts["G"]), words, sym))
assert np.isinf(A[0, 1])        # 0 -> 1: gains cannot resurrect a lost allele
assert A[1, 0] == 1 and A[1, 2] == 1
K = minplus(A.T, A)             # shared-ancestor kernel, zero WGDs
assert K[0, 1] == 1             # ancestor 1 loses one copy toward 0
print("PASS step 2")
```

## Step 3: model assets (copy, don't type)
**What/why.** `ancestral_parts.txt` holds the exported primitive transducers (`LOH`, `LG`, `L_LOH`, `G`, `A = L_LOH·G`,
compact `AD`). The three `generated/*.hpp` files are admissible lower-bound tables. Both engines read them.
```bash
REF=~/Projects/project-medicc/medicc-full-fst-lloh-wgd2-capped-x2/medicc/cext/full_fst
NEW=~/Projects/project-medicc/feature-full-fst-capped-wgd2/medicc/cext/full_fst
for e in ancestors direct_pair; do
  mkdir -p $NEW/$e/assets $NEW/$e/src/generated
  cp $REF/$e/assets/ancestral_parts.txt $NEW/$e/assets/
  cp $REF/$e/src/generated/*.hpp $NEW/$e/src/generated/
done
```
Read the first ~30 lines of `ancestral_parts.txt` to learn the format: a header `NAME n_states start n_arcs`, one line
of final weights, then `source input output weight target` rows. Symbol 9 is `X`.
**Test:** `tests/migration/test_03_parts_format.py`
```python
# <the header>
path = ROOT / "medicc/cext/full_fst/direct_pair/assets/ancestral_parts.txt"
lines = iter(path.read_text().splitlines())
blocks = {}
for header in lines:
    name, n, start, count = header.split()
    finals = list(map(int, next(lines).split()))
    arcs = [tuple(map(int, next(lines).split())) for _ in range(int(count))]
    assert len(finals) == int(n)
    blocks[name] = (int(n), int(start), arcs)
assert list(blocks) == ["LOH", "LG", "L_LOH", "G", "A", "AD"]
assert (blocks["A"][0], len(blocks["A"][2])) == (212, 15688)
assert (blocks["AD"][0], len(blocks["AD"][2])) == (74, 3108)
# L_LOH: state = number of active loss intervals; X resets to state 0 for free.
for s, a, b, w, t in blocks["L_LOH"][2]:
    if a == 9:
        assert (b, w, t) == (9, 0, 0)
    else:
        assert b == max(0, a - t) and w == max(0, t - s)
same = (ROOT / "medicc/cext/full_fst/ancestors/assets/ancestral_parts.txt").read_text()
assert same == path.read_text()
print("PASS step 3")
```

## Step 4: `scripts/generate_chromosome_bound.py`: regenerate the chromosome bound tables
**What/why.** Computes per-chromosome lower bounds for the capped model (`3*ka+kb` WGD branches). It reads the
*existing* header to keep its hand-written DP suffix, so Step 3 must come first.
**REF:** `scripts/generate_chromosome_bound.py`. **Type exactly, with one deviation:** in `main()` the engine list is
`['ancestors', 'direct_pair']` (no `reference_pair`). Change the final print to "both engines".
**Test (shell):**
```bash
python -B scripts/generate_chromosome_bound.py
for e in ancestors direct_pair; do
  cmp medicc/cext/full_fst/$e/src/generated/chromosome_bound.hpp \
      ~/Projects/project-medicc/medicc-full-fst-lloh-wgd2-capped-x2/medicc/cext/full_fst/$e/src/generated/chromosome_bound.hpp || exit 1
done && echo PASS step 4
```
(In REF this regeneration reproduces the committed file byte for byte. I verified that.)

## Step 5: `scripts/generate_capped_ad.py`: the compact AD factor
**What/why.** OpenFST composes `A` with capped `D`, then encodes, determinises and minimises it to 74 states. This
replaces `A D A D A` by `(AD)(AD)A`. It is a representation change, not an approximation.
**REF:** `scripts/generate_capped_ad.py`. **Type exactly.**
**Test (shell, about 1 minute):**
```bash
python -B scripts/generate_capped_ad.py --tool . \
   --source medicc/cext/full_fst/direct_pair/assets/ancestral_parts.txt --output /tmp/AD.txt
# expect JSON with "compact": {"arcs": 3108, "states": 74} and both equivalence flags true
sed -n '17609,$p' medicc/cext/full_fst/direct_pair/assets/ancestral_parts.txt | cmp - /tmp/AD.txt && echo PASS step 5
```

---

# Phase 2: `direct_pair` native engine (pair distances + event trace)
Needs only `g++`, no OpenFST or Boost. Read guide pages 5–7 and 11–12 first. All files go in
`medicc/cext/full_fst/…` in NEW. **Header check** (used below) for a header `F`:
```bash
g++ -std=c++17 -fsyntax-only -x c++ -I$CONDA_PREFIX/include F && echo OK
```
Some headers are *textually included in the middle* of another file and can't compile alone. Their check happens at the
step that includes them, as noted.

## Step 6: `common/src/wgd_model.hpp`: the capped doubling, in one function
**What/why.** `double_cn(c) = c ≤ 4 ? 2c : 8`. This one line *is* the "capped-x2" in the branch name.
**REF:** same path, 5 lines. **Type exactly.**
**Test:** `tests/migration/cpp/test_06_wgd_model.cpp`, run it with
`g++ -std=c++17 tests/migration/cpp/test_06_wgd_model.cpp -o /tmp/t06 && /tmp/t06 && echo PASS step 6`
```cpp
#include "../../../medicc/cext/full_fst/common/src/wgd_model.hpp"
static_assert(wgd_model::double_cn(0) == 0);
static_assert(wgd_model::double_cn(3) == 6);
static_assert(wgd_model::double_cn(4) == 8);
static_assert(wgd_model::double_cn(5) == 8);
static_assert(wgd_model::double_cn(8) == 8);
int main() { return 0; }
```

## Step 7: `common/src/feasible_path.hpp`: greedy feasible path → upper bound
**What/why.** Walks one legal path to get an *upper bound* U for the exact search (guide p.12). A failed walk proves
nothing. It is included after `Side` is defined, so it can't compile alone.
**REF:** same path (40 lines). **Type exactly.** **Test:** compiled in Step 16 (and in Step 28 for ancestors).

## Step 8: `direct_pair/src/budget_2m_core.hpp`: storage, primitive reader, DFA minimiser
**What/why.** Defines `Arc`, `Part` (a small primitive transducer), `Node` / `DFA` (the deterministic acceptor),
`readpart` (parses `ancestral_parts.txt`), `linear`, and the acyclic `minimize`. The "2m" in the name is historical;
**it is not a state cap**.
**REF:** `direct_pair/src/budget_2m_core.hpp` (292 lines). **Type exactly.** **Test:** header check → `OK`.

## Step 9: `direct_pair/src/shared_engine.hpp`: model lifetime and C ABI basics
**What/why.** `struct AM` (the loaded model: owns all `Part`s and the scalar score cache), `am_load` / `am_free` /
`am_string` / `am_import` / `am_size`… and the weighted `unite`. Read the ownership comment at ~line 246.
**REF:** `direct_pair/src/shared_engine.hpp` (484). **Type exactly.** **Test:** header check.

## Step 10: `direct_pair/src/shared_trace.hpp`: shared traceback helpers
**REF:** same (115). **Type exactly.** **Test:** header check.

## Step 11: `direct_pair/src/lazy_engine.hpp`: older lazy weighted-subset cascade
**What/why.** Not the main path any more, but `nfa_engine.hpp` uses types from `lazy_cascade::`.
**REF:** same (468). **Type exactly.** **Test:** header check.

## Step 12: `direct_pair/src/lazy_trace.hpp`
**REF:** same (118). **Type exactly.** **Test:** header check.

## Step 13: `direct_pair/src/nfa_engine.hpp`: primitive-tuple A* (scalar-score flavour)
**What/why.** This is the `direct_pair` copy of the join engine. It differs from the `ancestors` copy (Step 27): pair
scoring needs only the optimum, so it rejects transitions that merely *tie* the cutoff and does not rebuild the tied
language. Read guide p.11–13 while typing.
**REF:** `direct_pair/src/nfa_engine.hpp` (743). **Type exactly.** **Test:** header check.

## Step 14: `direct_pair/src/direct_trace.hpp`: one-sided A* for directed scores
**What/why.** Minimises `T(p→c) + M(c)` for a fixed parent. The heavily commented C++ here is the best introduction to
the code style (indices instead of pointers, borrowed references).
**REF:** same (199). **Type exactly.** **Test:** header check.

## Step 15: `direct_pair/src/event_trace.hpp`: recover intermediate profiles of one optimal history
**What/why.** Exports `am_pair_history`: JSON `{"cost", "wgds", "stages"}` on the same `(AD)^k A` cascade, preferring
fewer WGDs. It is included at the *end* of `pairwise.cpp`.
**REF:** same (75). **Type exactly.** **Test:** compiled in Step 16.

## Step 16: `direct_pair/src/pairwise.cpp`: exact factorised pair scorer (the library entry)
**What/why.** `am_pair_score(h, x, y, directed, stats)`. It picks global WGD counts `(ka, kb)` first, then scores
`X`-blocks short-first with witness upper bounds (`witness`, `chromosome_witness`) and caching. Guide p.7 explains why
the min over `k,l` must stay *outside* the per-chromosome sum.
**REF:** same (234). **Type exactly.** Then build it by hand once:
```bash
D=medicc/cext/full_fst/direct_pair; mkdir -p $D/build && cp $D/assets/ancestral_parts.txt $D/build/
g++ -O3 -std=c++17 -pthread -shared -fPIC -DANCESTRAL_MAX_WGDS=2 $D/src/pairwise.cpp -o $D/build/libpair_direct.so
```
**Test:** `tests/migration/test_16_raw_ctypes.py`. This is the raw C ABI with no wrapper, so you see what
`_bindings.py` will hide.
```python
# <the header>
import ctypes as C
import math
B = ROOT / "medicc/cext/full_fst/direct_pair/build"
lib = C.CDLL(str(B / "libpair_direct.so"))
lib.am_load.argtypes = [C.c_char_p]; lib.am_load.restype = C.c_void_p   # without restype: 32-bit int, crash
lib.am_free.argtypes = [C.c_void_p]; lib.am_free.restype = None
lib.am_pair_score.argtypes = [C.c_void_p, C.c_char_p, C.c_char_p, C.c_int, C.POINTER(C.c_double)]
lib.am_pair_score.restype = C.c_double
lib.am_pair_history.argtypes = [C.c_void_p, C.c_char_p, C.c_char_p, C.c_int]
lib.am_pair_history.restype = C.c_char_p
h = lib.am_load(str(B / "ancestral_parts.txt").encode())
stats = (C.c_double * 6)()
s = lambda x, y, d=0: lib.am_pair_score(h, x.encode(), y.encode(), d, stats)
assert s("1", "2") == 1 and s("11", "22") == 1
assert s("0", "1") == 1 and math.isinf(s("0", "1", 1))   # kernel vs directed
assert s("22", "11", 1) == 1 and s("2X2", "1X1", 1) == 2
assert s("4X5", "8X8", 1) == 1                           # one capped WGD
assert s("1X1", "8X8", 1) == 4                           # 2 gains + 2 WGDs
assert s("00014", "13207") == 5 and s("0X0", "7X7") == 8
assert lib.am_pair_history(h, b"1X1", b"8X8", 4) == b'{"cost":4,"wgds":2,"stages":["4X4","8X8","8X8"]}'
lib.am_free(h)
print("PASS step 16")
```

---

# Phase 3: build script and pair-scorer Python API

## Step 17: `scripts/build.py`: build the native engines
**What/why.** Compiles each engine into `medicc/cext/full_fst/<engine>/build/` and copies the asset beside it.
**REF:** `scripts/build.py`. **Type it with these deviations:**
* engine list: `ancestors` and `direct_pair` only. **Comment out the `ancestors` tuple until Step 30.**
* drop `--legacy` and `build_legacy.py` (master builds fstlib with `python setup.py build_ext --inplace`), and drop the
  "Missing Python extensions" check.
* drop `--record-only` and the `code_manifest.json` / `build_manifest.json` hashing (`provenance.py` is not migrated).

**Test:**
```bash
rm -rf medicc/cext/full_fst/direct_pair/build && python -B scripts/build.py && python -B tests/migration/test_16_raw_ctypes.py
```
(`build/` is already ignored by master's `.gitignore`.)

## Step 18: `medicc/full_fst/_bindings.py`: all ctypes signatures in one place
**What/why.** Library paths, `argtypes` / `restype` for every exported function, and `_check_model` (refuses a library
that isn't the L_LOH / 2-WGD / fused build).
**REF:** same. **Type exactly.** Also create an **empty** `medicc/full_fst/__init__.py` (filled in Step 46).
**Test:** `tests/migration/test_18_bindings.py`
```python
# <the header>
from medicc.full_fst import _bindings as b
lib = b.load_pair_library("direct_pair")
assert (lib.am_max_wgds(), lib.am_lloh_model(), lib.am_fused_model()) == (2, 1, 1)
assert lib.am_pair_score_only() == 1 and lib.am_pair_search_uncapped() == 1
assert b.build_directory("direct_pair").name == "build"
try:
    b.load_pair_library("reference_pair")      # not migrated -> no library file
except OSError:
    pass
else:
    raise AssertionError("reference_pair unexpectedly present")
print("PASS step 18")
```

## Step 19: `medicc/full_fst/pair_scorer.py`: `PairScorer` (exact pair / directed scores)
**What/why.** Owns one native model and its cache, one per thread. `score(x, y)` = kernel `K`,
`score(x, y, directed=True)` = `T`. `clear()` resets the cache at row boundaries.
**REF:** same. **Type exactly** (the `reference_pair` branch of the backend check may stay; it just fails to load).
**Test:** `tests/migration/test_19_pair_scorer.py`
```python
# <the header>
import math
from medicc.full_fst.pair_scorer import Model as PairScorer

with PairScorer() as m:
    assert m.backend == "direct_pair"
    assert m.score("1X1", "8X8") == m.score("8X8", "1X1")            # K is symmetric
    assert m.score("1X1", "8X8", directed=True) == 4
    assert math.isinf(m.score("8X8", "1X1", directed=True)) is False  # losses are fine
    assert math.isinf(m.score("0", "1", directed=True))
    assert m.score("12", "1X2") == math.inf                          # grid mismatch
    m.score("1122", "2211"); assert len(m.stats()) == 6
    m.clear(); assert m.stats() == [0.0] * 6
    for bad in ["9", "1a"]:
        try: m.score(bad, bad)
        except ValueError: pass
        else: raise AssertionError(bad)
m.close()                                                             # idempotent
try: m.score("1", "1")
except ValueError: print("PASS step 19")
```

## Step 20: `tests/check_pairs.py`: the exhaustive oracle comparison
**What/why.** Compares every pair of 2-bin words (`##` and `#X#`) against the oracle: 26,244 comparisons on one backend.
**REF:** `tests/check_pairs.py`. **Type exactly, except** `models = [Model(backend="direct_pair")]` and the PASS text.
**Test:** `python -B tests/check_pairs.py` → `PASS: 26244 oracle comparisons…` (REF prints 52,488 because it has two backends).

---

# Phase 4: `ancestors` native engine (all-tie join + resolve; needs OpenFST + Boost)
Read guide pages 8–16 before starting. The header check now also needs Boost (it's in `$CONDA_PREFIX/include`).

## Step 21: `common/src/automata_types.hpp`: shared `Node` / `DFA` / `INF`
**What/why.** In `direct_pair` these types live inside `budget_2m_core.hpp`. The ancestors engine shares them with the
OpenFST adapter, so they get their own header.
**REF:** same (21). **Type exactly.**
**Test:** `tests/migration/cpp/test_21_types.cpp` (`g++ -std=c++17 … -o /tmp/t21 && /tmp/t21`)
```cpp
#include "../../../medicc/cext/full_fst/common/src/automata_types.hpp"
#include <cassert>
int main() {
    Node n; DFA d;
    for (int c = 0; c < 10; c++) assert(n.to[c] == -1 && n.w[c] == 0);
    assert(n.fin == INF && d.start == -1 && d.offset == 0 && d.v.empty());
    return 0;
}
```

## Step 22: `common/src/openfst_language.hpp` + `.cpp`: the only OpenFST translation unit
**What/why.** Converts MEDICC DFAs ↔ OpenFST acceptors (labels shifted by +1, because 0 is epsilon). Provides
determinize, weighted union and rmepsilon for the tied ancestor languages. Kept separate so the search kernels never see
OpenFST templates.
**REF:** both files (8 + 106). **Type exactly.**
**Test:**
```bash
g++ -std=c++17 -fPIC -c -I$CONDA_PREFIX/include medicc/cext/full_fst/common/src/openfst_language.cpp -o /tmp/ofl.o && echo PASS step 22
```

## Step 23: `common/src/boost_search.hpp`: Boost.Graph A* adapter
**What/why.** Adapts the implicit, lazily generated search graph to `astar_search_no_init`, with an exception used to
stop at a proven bound and stable tie ordering (guide p.13, p.21).
**REF:** same (122). **Type exactly.** **Test:** header check → `OK`.

## Step 24: `ancestors/src/budget_2m_core.hpp`
**What/why.** Same role as Step 8. It differs from the `direct_pair` copy in 20 lines: it includes
`automata_types.hpp` / `openfst_language.hpp` instead of defining `Node` / `DFA` / `INF` itself.
**REF:** `ancestors/src/budget_2m_core.hpp`. Type it fresh, or `cp` the direct_pair file and type only the differences
(`diff REF/direct_pair/src/budget_2m_core.hpp REF/ancestors/src/budget_2m_core.hpp`).
**Test:** header check, then `cmp` against REF.

## Step 25: `ancestors/src/shared_engine.hpp`
**What/why.** Same as Step 9, but `unite` calls OpenFST (`library_language::united`) instead of the hand-written subset
construction, and there is no scalar cache.
**REF:** same (406). **Type** (or cp + diff as in Step 24). **Test:** header check + `cmp`.

## Step 26: identical support headers (copy)
`shared_trace.hpp`, `lazy_engine.hpp`, `lazy_trace.hpp` are byte-identical to your `direct_pair` copies:
```bash
for f in shared_trace lazy_engine lazy_trace; do
  cp medicc/cext/full_fst/direct_pair/src/$f.hpp medicc/cext/full_fst/ancestors/src/
  cmp medicc/cext/full_fst/ancestors/src/$f.hpp ~/Projects/project-medicc/medicc-full-fst-lloh-wgd2-capped-x2/medicc/cext/full_fst/ancestors/src/$f.hpp
done && echo PASS step 26
```

## Step 27: `ancestors/src/pattern_bound.hpp`: sibling pattern lower bound
**What/why.** Projects each primitive state onto 4 classes (free decrease / increase / both / neither) and solves the
relaxed sibling problem. This gives a safe lower bound (guide p.15). It is included in the middle of `nfa_engine.hpp`
(line 227).
**REF:** same (162). **Type exactly.** **Test:** compiled in Step 29.

## Step 28: `ancestors/src/segmented_join.hpp`: exact fallback at safe `X` cuts
**What/why.** After 250k expansions, restart as independent segments when both child DFAs have a single state at an
`X` layer (guide p.16). This selects an algorithm; it is **not** a limit. It is included at `nfa_engine.hpp:510`.
**REF:** same (131). **Type exactly.** **Test:** compiled in Step 29.

## Step 29: `ancestors/src/nfa_engine.hpp`: the all-tie bottom-up join (the heart of the branch)
**What/why.** `am_nfa_join`: forward A* → best cost; backward pass keeps **every tight arc**; then determinize, unite
WGD languages and minimise. Also defines `ExecutionBudget` (`am_budget_*`) and `am_choose_profile`. While typing,
compare with your Step 13 copy: each "score-only" shortcut there is deliberately absent here.
**REF:** same (761). **Type exactly.**
**Test:** `g++ -std=c++17 -fsyntax-only -x c++ -I$CONDA_PREFIX/include medicc/cext/full_fst/ancestors/src/nfa_engine.hpp && echo OK`.
This also compiles Steps 7, 27 and 28.

## Step 30: `ancestors/src/ancestors.cpp`: top-down resolve + library entry
**What/why.** `am_direct_resolve(h, parent, message, upper, &cost, stats)`: one-sided A* choosing one optimal child
profile. It deliberately keeps the original `A, D, A` order, so tie choices stay historical (guide p.14).
**REF:** same (179). **Type exactly.** Then **uncomment the `ancestors` tuple in `scripts/build.py`** and run
`python -B scripts/build.py`.
**Test:** `tests/migration/test_30_ancestor_abi.py`
```python
# <the header>
from medicc.full_fst._bindings import load_ancestor_library, build_directory
lib = load_ancestor_library()                         # also runs _check_model
h = lib.am_load(str(build_directory("ancestors") / "ancestral_parts.txt").encode())
assert h, lib.am_error()
d = lib.am_string(b"12X3")
assert lib.am_size(d) == 5 and lib.am_offset(d) == 0  # linear acceptor: len+1 states
assert lib.am_start(d) >= 0
lib.am_drop(d); lib.am_free(h)
print("PASS step 30")
```

---

# Phase 5: messages and the `Model` API

## Step 31: `medicc/full_fst/message.py`: Python owner of one native DFA
**What/why.** `Message` owns the pointer (`close()` calls `am_drop` exactly once) and keeps its `Model` alive.
`language()` enumerates words (tiny tests only), `to_fst()` exports to fstlib, and `choose()` gives one cheapest word.
**REF:** same. **Type exactly.**
**Test:** `tests/migration/test_31_message.py`. It uses a stand-in object with a `.lib` attribute, so you see the
ownership contract without `Model`.
```python
# <the header>
from types import SimpleNamespace
from medicc.full_fst._bindings import load_ancestor_library
from medicc.full_fst.message import Message
owner = SimpleNamespace(lib=load_ancestor_library())
m = Message(owner, owner.lib.am_string(b"12"))
assert m.language() == {"12": 0} and m.size == 3 and m.offset == 0
assert m.choose() == ("12", 0)
m.close(); m.close()                                  # second close is a no-op
try: m.size
except ValueError: pass
else: raise AssertionError("closed message usable")
try: Message(owner, None)                             # NULL handle -> native error text
except RuntimeError: print("PASS step 31")
```

## Step 32: `medicc/full_fst/model.py`: `Model.string` / `join` / `resolve`, plus `ExecutionBudget`
**What/why.** The three operations from guide p.6. `join` computes an exact upper bound from the children's chosen
words before calling `am_nfa_join`. `resolve` does the same with a directed score. `from_fst` imports a fixed-length
deterministic acceptor. One `Model` per thread.
**REF:** same. **Type exactly.**
**Test:** `tests/migration/test_32_model.py` (all values checked against REF)
```python
# <the header>
import medicc
from medicc.full_fst.model import Model, ExecutionBudget

m = Model()
one, two = m.string("1"), m.string("2")
j = m.join(one, two)
assert j.language() == {"1": 1, "2": 1} and j.offset == 1   # BOTH ancestors survive (guide p.10)
assert m.resolve("1", j) == ("1", 1)
assert m.resolve("2", j) == ("2", 1)
assert m.resolve("3", j) == ("2", 2)                         # cost = edge 1 + stored message cost 1
j2 = m.join(m.string("11"), m.string("22"))
assert j2.language() == {"11": 1, "22": 1}                   # NOT {"12","21"}: correlations kept
sym = medicc.create_symbol_table(8, "X")
back = m.from_fst(j2.to_fst(sym), 2)                         # fstlib round trip
assert back.language() == j2.language()
m.execution_budget = ExecutionBudget(2)
assert m.join(one, two).language() == {"1": 1, "2": 1}
assert {"nfa_expansions", "down_expansions"} <= set(m.stats())
for bad in ["9", "1a"]:
    try: m.string(bad)
    except ValueError: pass
    else: raise AssertionError(bad)
try: Model(max_wgds=1)
except ValueError: pass
else: raise AssertionError("max_wgds=1 accepted")
m.close(); m.close()
print("PASS step 32")
```

## Step 33: port two native-language checks from REF
**What/why.** `check_resolve_tie.py` pins the historical tie choice (`5448 → 2264`, cost 8).
`check_ancestor_bounds.py` compares *complete weighted languages* of joins against the oracle (guide p.22).
**REF:** `tests/check_resolve_tie.py`, `tests/check_ancestor_bounds.py`. **Type exactly** into NEW `tests/`.
**Test:** `python -B tests/check_resolve_tie.py` → `{"status": "ok", …}`, then `python -B tests/check_ancestor_bounds.py`
→ PASS line.

---

# Phase 6: tree workflow

## Step 34: `medicc/full_fst/shortening.py`: collapse equal adjacent column pairs (pairs only)
**What/why.** Precomputes per-sample change masks. A pair keeps the union of both masks and every `X`. Used only for
distance scoring; messages keep all bins.
**REF:** same (22). **Type exactly.**
**Test:** `tests/migration/test_34_shortening.py`
```python
# <the header>
from medicc.full_fst import shortening
from medicc.full_fst.pair_scorer import Model
a, b = shortening.prepare("1122X33"), shortening.prepare("1222X33")
assert shortening.pair(a, b) == ("112X3", "122X3")
with Model() as m:
    assert m.score(*shortening.pair(a, b)) == m.score("1122X33", "1222X33")
try: shortening.pair(a, shortening.prepare("1"))
except ValueError: print("PASS step 34")
```

## Step 35: `medicc/full_fst/scheduler.py`: bottom-up joins, then top-down resolves, on threads
**What/why.** Input is plain dicts (`children`, `leaves`, `root`, `normal`). The ready-queue schedule starts a node as
soon as its children are done. Each worker thread gets its own `Model`. It returns `(profiles, messages, metrics)`.
**REF:** same. **Type exactly.**
**Test:** `tests/migration/test_35_scheduler.py`
```python
# <the header>
from medicc.full_fst.scheduler import reconstruct
leaves = {"a": "11X1", "b": "22X2", "c": "21X1"}
children = {"r": ["a", "s"], "s": ["b", "c"]}
expected = {"a": "11X1", "b": "22X2", "c": "21X1", "r": "11X1", "s": "11X1"}
for scheduling in ["ready", "levels"]:
    profiles, messages, metrics = reconstruct(children, leaves, "r", "11X1", workers=2, scheduling=scheduling)
    assert profiles == expected, (scheduling, profiles)
    assert messages == {} and metrics["up_nodes"] == metrics["down_nodes"] == 2
for bad in [dict(children={"r": ["a", "b", "c"]}), dict(children={"r": ["a", "zz"]})]:
    try: reconstruct(bad["children"], leaves, "r", "11X1")
    except ValueError: pass
    else: raise AssertionError(bad)
print("PASS step 35")
```

## Step 36: `medicc/full_fst/core.py`: encode → distances → tree → ancestors → branch lengths
**What/why.** The public workflow (`Result`, `encode`, `distance_matrix`, `normalize_tree`, `infer_tree`, `decode`,
`reconstruct`).
**REF:** same. **Type it with these deviations (B2):**
1. `encode`: keep master's sample order. Replace the natural-sort `labels = [normal_name] + sorted(...)` with
   `labels = list(input_df.index.get_level_values("sample_id").unique())`. Keep every check.
2. `normalize_tree`: **master's naming scheme** (see box below). Change REF's renaming loop so that:
   * it walks `root.find_clades(order="postorder")` (deepest first, like NJ's join order) instead of
     `get_nonterminals()` (preorder);
   * the new name is `f"internal{counter}"` (1-based, no underscore), skipping names already in `used`;
   * existing names are never touched unless `rename=True`. The only call with `rename=True` is FastME's in
     `infer_tree`, because FastME's internal labels are meaningless.
3. `reconstruct(...)`: add `tree_method="nj"`. When `input_tree is None`:
   `"nj"` → `normalize_tree(medicc.core.infer_tree_topology(distances.values, distances.index, normal_name), normal_name, labels)`
   (master's NJ, **no `rename`**, so its `internalN` names survive unchanged). `"fastme"` → REF's `infer_tree(...)`
   call. Anything else → `ValueError`. Record `details["tree_method"]`.

> **Internal node names = master's scheme.** master names come from `medicc/nj.py`: `internal1`, `internal2`, … in
> NJ join order (1-based, no underscore). A supplied `--tree` keeps the user's names. NEW must give:
> * NJ: names identical to legacy master. Verified: on `simple_example`, full model + master NJ gives *exactly*
>   `internal1={taxon_2,taxon_3}`, `internal4={+taxon_5}`, `internal3={+taxon_1}`, `internal2={all tumours}`, the same
>   as legacy `medicc2`.
> * supplied tree: the user's names; only unnamed nodes get `internalN`.
> * FastME: `internalN`, numbered bottom-up.
> * NNI keeps every name (moves only swap subtrees, Step 41). The output wrapper adds an *unnamed* root, as in master.
> * never `inferred_N` or `internal_N`. (master's `internal_1` with an underscore only appears in its
>   one-tumour special case, which the full model rejects: it needs at least two tumours.)
4. The `if nni:` block: **leave a `raise NotImplementedError("Step 43")` for now.** You'll write it in Step 43.

Keep `infer_tree` exactly as REF (FastME, long-name remap).
**Test:** `tests/migration/test_36_core.py`
```python
# <the header>
import tempfile, os
from io import StringIO
import pandas as pd
from Bio import Phylo
import medicc
from medicc.full_fst.core import encode, decode, distance_matrix, reconstruct

rows = []
cn_a = {"a": ("22", "11"), "b": ("22", "12"), "c": ("11", "22"), "d": ("21", "22")}
for name, (chr1, chr2) in cn_a.items():
    for chrom, word in [("chr1", chr1), ("chr2", chr2)]:
        for i, c in enumerate(word):
            rows.append((name, chrom, i * 10, i * 10 + 9, int(c), 1))
tmp = tempfile.mkdtemp()
path = os.path.join(tmp, "toy.tsv")
pd.DataFrame(rows, columns=["sample_id", "chrom", "start", "end", "cn_a", "cn_b"]).to_csv(path, sep="\t", index=False)
frame = medicc.io.read_and_parse_input_data(path, normal_name="normal")

genomes, coords, cols = encode(frame, "normal")
assert genomes["a"] == "22X11X11X11" and genomes["normal"] == "11X11X11X11"   # allele-major blocks
assert decode(genomes, coords, cols).astype(str).equals(frame.astype(str).sort_index())
d, stats = distance_matrix(genomes, workers=2)
assert (d.values == d.values.T).all() and d.loc["a", "b"] == 1 and d.loc["a", "c"] == 2
tree = Phylo.read(StringIO("((a,b),(c,d));"), "newick")
r = reconstruct(frame, normal_name="normal", input_tree=tree, n_cores=2)
assert r.details["tree_length"] == 4 and r.pairwise_distances is None
assert r.tree.root.name == "normal" and len(r.tree.root.clades) == 1
clades = lambda t: {n.name: sorted(x.name for x in n.get_terminals())
                     for n in t.get_nonterminals() if n.name != "normal"}
assert clades(r.tree) == {"internal1": ["a", "b"], "internal2": ["c", "d"], "internal3": ["a", "b", "c", "d"]}  # bottom-up (REF: preorder inferred_1..3)
named = Phylo.read(StringIO("((a,b)left,(c,d));"), "newick")                     # user names are kept
assert set(clades(reconstruct(frame, normal_name="normal", input_tree=named).tree)) == {"left", "internal1", "internal2"}
r = reconstruct(frame, normal_name="normal", n_cores=2)                          # NJ default
assert r.details["tree_method"] == "nj" and r.details["tree_length"] == 4
assert clades(r.tree) == {"internal1": ["a", "b"], "internal3": ["c", "d"], "internal2": ["a", "b", "c", "d"]}  # master NJ names
r = reconstruct(frame, normal_name="normal", n_cores=2, tree_method="fastme",
                fastme=os.path.join(os.environ["CONDA_PREFIX"], "bin", "fastme"))
assert r.details["tree_length"] == 4
assert sorted(clades(r.tree)) == ["internal1", "internal2", "internal3"]
print("PASS step 36")
```

---

# Phase 7: exact event calling

## Step 37: `medicc/full_fst/events.py`: one replayable minimum-cost history per branch
**What/why.** The native trace gives WGD count + stage profiles. A small OpenFST shortest path splits each stage into
L_LOH / G active-interval counts, which become unit events. `replay` must reach the child, and the event count must
equal the score (guide p.19). `call_tree_events` builds a master-shaped event table plus extra columns
(`parent_id`, `event_order`, `stage`, `cost`, `genome_*`).
**REF:** same. **Type exactly.** In `core.reconstruct`, the `if events:` block is already typed from Step 36.
**Test:** `tests/migration/test_37_events.py`
```python
# <the header>
from medicc.full_fst.events import EventCaller, Event, replay
assert replay("2X2", [Event("loss", 0, 1, 0), Event("wgd", -1, -1, 0)]) == "2X4"
with EventCaller() as c:
    h = c.branch("1X1", "8X8")
    assert (h.cost, h.wgds) == (4, 2)
    assert [(e.type, e.start, e.end, e.stage) for e in h.events] == [
        ("gain", 0, 1, 0), ("gain", 2, 3, 0), ("wgd", -1, -1, 0), ("wgd", -1, -1, 1)]
    assert [(e.type, e.start, e.end) for e in c.branch("22", "11").events] == [("loss", 0, 2)]
    assert [e.type for e in c.branch("4X5", "8X8").events] == ["wgd"]
    assert c.branch("11", "11").events == []
    for x, y in [("1X1", "111"), ("9", "8"), ("", ""), ("0", "1")]:
        try: c.branch(x, y)
        except ValueError: pass
        else: raise AssertionError((x, y))
print("PASS step 37 (part 1)")
```
Then add (re-using the toy frame from Step 36) `reconstruct(frame, normal_name="normal", input_tree=tree, events=True)`
and assert `len(r.events) == 4`, `set(r.events["type"]) == {"gain"}`, and that the rows have
`parent_id` starting with `"internal"` for the two clade gains.

## Step 38: `medicc/event_reconstruction.py` (edit): keep every WGD in the overlaps (bug fix, see D3)
**What/why.** master overwrites one row label for each WGD, so only the last WGD survives in `_events_overlap.tsv`.
With two WGDs allowed per edge, this now happens easily.
**Order:** type the test **first**, run it on the unmodified file and watch it fail (this proves the bug), then port the fix.
**REF:** `medicc/event_reconstruction.py:399-407`. The `wgd_rows` list + `pd.concat` block replaces master's
`for ind, _ in …iterrows(): all_events.loc[('all','0','0')] = …` loop. Port only that block.
(The `cn_child … .astype(float)` line is part of Phase 0's toolchain changes; the reworded messages are not ported.)
**Test:** `tests/migration/test_38_overlap_wgd.py` (needs pyranges, which is in both envs)
```python
# <the header>
import pandas as pd
from medicc.event_reconstruction import overlap_events

rows = [("branch_a", "chr0", 0, 99, "both", "wgd", 0),
        ("branch_b", "chr0", 0, 99, "both", "wgd", 0),
        ("branch_b", "chr1", 0, 99, "cn_a", "gain", 2)]
events = pd.DataFrame(rows, columns=["sample_id", "chrom", "start", "end", "allele", "type", "cn_child"]
                      ).set_index(["sample_id", "chrom", "start", "end"])
out = overlap_events(events_df=events, alleles=["cn_a", "cn_b"], normal_name="normal",
                     chromosome_bed="default", regions_bed="none").reset_index()
wgd = out[out["event"] == "WGD"]
assert sorted(wgd["branch"]) == ["branch_a", "branch_b"], wgd   # master: only ['branch_b']
print("PASS step 38")
```
(Checked: master's code returns one WGD row, `branch_b`; REF's code returns both.)
---

# Phase 8: NNI

## Step 39: `medicc/tools.py` (edit): add `sum_of_branch_length`
**REF:** last 4 lines of `medicc/tools.py`. Add only this function. **Keep master's logger name.**
**Test:** one-liner in `tests/migration/test_39_tools.py`: on `Phylo.read(StringIO("((a:1,b:2):3,c):0;"), "newick")`,
assert `medicc.tools.sum_of_branch_length(t) == 6`.

## Step 40: `medicc/tree_hash.py`: rotation-invariant topology hash
**REF:** same (from your pre-codex work). **Type exactly.**
**Test:** `tests/migration/test_40_tree_hash.py`
```python
# <the header>
from io import StringIO
from Bio import Phylo
from medicc.tree_hash import get_topology_hash, get_canonical_newick
t = lambda s: Phylo.read(StringIO(s), "newick")
assert get_topology_hash(t("((a,b),(c,d));")) == get_topology_hash(t("((d,c),(b,a));"))
assert get_topology_hash(t("((a,b),(c,d));")) != get_topology_hash(t("((a,c),(b,d));"))
assert get_canonical_newick(t("((b,a)x,(d,c)y);")) == get_canonical_newick(t("((c,d),(a,b));"))
print("PASS step 40")
```

## Step 41: `medicc/nni.py`: topology-only NNI moves
**What/why.** `_NNIMove`, `_apply_nni_swap`, `_enumerate_moves`, `nni_neighbors`: two moves per internal edge
(u, v), swapping the uncle with one of v's children.
**REF:** `medicc/nni.py` lines **1–89 only**. Do **not** type `_eval_nni_neighbor` / `evaluate_nni_neighbors_parallel`
here; they move to Step 43 without FST arguments.
**Test:** `tests/migration/test_41_nni_moves.py`
```python
# <the header>
from io import StringIO
from Bio import Phylo
from medicc.full_fst.core import normalize_tree
from medicc.nni import _enumerate_moves, nni_neighbors
from medicc.tree_hash import get_topology_hash
labels = ["normal", "a", "b", "c", "d"]
t = normalize_tree(Phylo.read(StringIO("((a,b),(c,d));"), "newick"), "normal", labels, rename=True)
assert t.root.name == "normal" and len(_enumerate_moves(t)) == 4
t3 = normalize_tree(Phylo.read(StringIO("(a,(b,c));"), "newick"), "normal", ["normal", "a", "b", "c"], rename=True)
assert len(_enumerate_moves(t3)) == 2
hashes = {get_topology_hash(n) for n, _ in nni_neighbors(t)}
assert get_topology_hash(t) not in hashes                          # every move changes the topology
names = sorted(c.name for c in t.find_clades())
assert all(sorted(c.name for c in n.find_clades()) == names for n, _ in nni_neighbors(t))   # names preserved
print("PASS step 41")
```

## Step 42: `medicc/full_fst/ancestors.py`: tree reconstruction with message cache (no FST arguments)
**What/why.** Bottom-up joins into a `Cache` of `CachedMessage`s, then top-down resolves. For an NNI neighbour, only
the path from `nni_move.v_name` to the root is *dirty*. Resolves are reused only when both the message and the parent
profile are identical.
**REF:** `medicc/full_fst/ancestors.py`. **Type with these deviations (B1):**
* delete `FULL_TWO = object()` and the imports `sys`, `Path`, `medicc`, `medicc.ancestors`, `medicc.core`.
* `reconstruct(tree, samples_dict, normal_name, prune_weight=0, n_cores=None, upper_cache=False, old_uppass_cache=None, nni_move=None, **unused)`.
  The first check becomes only `if prune_weight != 0: raise ValueError("Requires prune(0)")`.
* `update_lengths(tree, ancestor_fsa, normal_name="diploid", log_progress=True)`: no `fst` argument and no FULL_TWO
  check. The default `"diploid"` matches master's `update_branch_lengths`.

**Test:** `tests/migration/test_42_incremental.py`
```python
# <the header>
import copy
from io import StringIO
from Bio import Phylo
from medicc.full_fst.core import normalize_tree
from medicc.full_fst.ancestors import Genome, reconstruct, update_lengths
from medicc.nni import nni_neighbors
from medicc.tools import sum_of_branch_length

s = {k: Genome(v) for k, v in {"normal": "11X11", "a": "22X11", "b": "22X12", "c": "11X22", "d": "21X22"}.items()}
t = normalize_tree(Phylo.read(StringIO("((a,b),(c,d));"), "newick"), "normal", list(s), rename=True)
anc, cache = reconstruct(t, s, "normal", upper_cache=True)
update_lengths(t, anc, "normal")
assert sum_of_branch_length(t) == 4
for neighbour, move in nni_neighbors(t):
    inc, _ = reconstruct(neighbour, s, "normal", old_uppass_cache=cache, nni_move=move)
    full = reconstruct(copy.deepcopy(neighbour), s, "normal")
    assert dict(inc) == dict(full), move.v_name                 # incremental == from scratch
for n_cores in (1, 2):
    assert dict(reconstruct(t, s, "normal", n_cores=n_cores)) == dict(anc)
try: reconstruct(t, s, "normal", prune_weight=1)
except ValueError: print("PASS step 42")
```

## Step 43: `medicc/full_fst/nni_search.py` (new file): neighbour evaluation + `nni_mode`, no FST arguments
**What/why.** Steepest-ascent NNI with a frontier, plateau handling, a visited-topology cache, and thread-parallel
neighbour evaluation over one shared, read-only message cache.
**REF sources** (type them into this one file):
* `medicc/nni.py:91-158`: `_eval_nni_neighbor`, `evaluate_nni_neighbors_parallel`
* `medicc/core.py:623-993`: `nni_mode`

**Deviations (B1):**
* module header: `import logging; logger = logging.getLogger(__name__)` (lands under master's `medicc` logger
  hierarchy), plus `from medicc.nni import _enumerate_moves, nni_neighbors, _apply_nni_swap`,
  `from medicc.tree_hash import get_topology_hash, get_canonical_newick`, `from medicc.tools import sum_of_branch_length`,
  `from .ancestors import reconstruct, update_lengths, validate_profiles`.
* `_eval_nni_neighbor(tree, move, old_uppass_cache, samples_dict, normal_name, prune_weight, visited_solutions, step=0)`
  and `evaluate_nni_neighbors_parallel(tree, old_uppass_cache, samples_dict, normal_name, prune_weight, n_cores, visited_solutions, step_start=0)`:
  no FST arguments. Inside, call `reconstruct(new_tree, samples_dict, normal_name, prune_weight=…, old_uppass_cache=…, nni_move=move, validate_result=False)`
  and `update_lengths(new_tree, ancestors, normal_name, log_progress=False)`.
* `nni_mode(tree, samples_dict, normal_name="diploid", prune_weight=0, nni_max_iter=20000, n_cores=None, nni_trace_dir=None, nni_max_frontier=10, nni_stop_on_plateau=False)`:
  no `upper_pass_fst` / `lower_pass_fst` / `la_nni`.
  * Delete both inner `_…_wrapper` functions. Call `reconstruct(..., upper_cache=True)` and
    `reconstruct(..., old_uppass_cache=..., nni_move=...)` directly.
  * Every `if la_nni: X else: Y` pair becomes one call.
  * `update_branch_lengths(t, fst, anc, normal)` → `update_lengths(t, anc, normal_name)`.
  * `medicc.ancestors.validate_ancestor_lengths(sel, normal)` → `validate_profiles(sel, normal_name)`.
  * `medicc.nni.X` → the imported names, and `medicc.tools.sum_of_branch_length` → `sum_of_branch_length`.

Then finish **`core.reconstruct`'s `if nni:` block** from REF `full_fst/core.py`. Replace the
`medicc.core.nni_mode(tree, genomes, FULL_TWO, FULL_TWO, …)` call with `nni_search.nni_mode(tree, {…Genome}, normal_name=…, …)`,
then use `resolve_tree(tree, genomes, normal_name, n_cores=workers)` and `update_lengths(tree, profiles, normal_name)`.
**Test:** `tests/migration/test_43_nni.py` (REF gives `[5.0, 4.0]` for both serial and parallel)
```python
# <the header>
import os, tempfile
from io import StringIO
from Bio import Phylo
from medicc.full_fst.core import normalize_tree
from medicc.full_fst.ancestors import Genome
from medicc.full_fst.nni_search import nni_mode
from medicc.tree_hash import get_topology_hash

s = {k: Genome(v) for k, v in {"normal": "11X11", "a": "22X11", "b": "22X12", "c": "11X22", "d": "21X22"}.items()}
start = normalize_tree(Phylo.read(StringIO("(a,(b,(c,d)));"), "newick"), "normal", list(s), rename=True)
target = get_topology_hash(Phylo.read(StringIO("(((a,b),(c,d)))normal;"), "newick"))
for n_cores in (1, 2):
    r = nni_mode(start, s, normal_name="normal", n_cores=n_cores, nni_max_iter=float("inf"),
                 nni_max_frontier=1, nni_stop_on_plateau=True, nni_trace_dir=tempfile.mkdtemp())
    assert r["trace"] == [5.0, 4.0], (n_cores, r["trace"])
    assert r["best_score"] == 4 and get_topology_hash(r["best_trees"][0]) == target
    assert r["step_records"] and r["step_records"][0][0] == 0
# ((a,c),(b,d)) is 2 rooted-NNI moves away: greedy descent stops at 5 (same as REF)
stuck = normalize_tree(Phylo.read(StringIO("((a,c),(b,d));"), "newick"), "normal", list(s), rename=True)
assert nni_mode(stuck, s, normal_name="normal", n_cores=1, nni_max_frontier=1, nni_stop_on_plateau=True)["trace"] == [5.0]
print("PASS step 43")
```
Also re-run Step 36's test with `reconstruct(frame, normal_name="normal", input_tree=Phylo.read(StringIO("(a,(b,(c,d)));"), "newick"), nni=True)`
→ `details["score_trace"] == [5.0, 4.0]`.

---

# Phase 9: CLI and packaging

## Step 44: `medicc2` → `medicc2-legacy` (move master's CLI out of the way, unchanged)
**What/why.** After this update the full model is MEDICC's default `medicc2`. master's OpenFST CLI stays available,
**byte for byte**, as `medicc2-legacy`. Everything in master that calls the script by name expects legacy behaviour, so
those callers move with it. Do this step *before* typing the new CLI, while `medicc2` is still master's script.

1. **Baseline first** (master's script, still called `medicc2`). Two identical runs give identical output, so `diff -r` works:
   ```bash
   python medicc2 examples/simple_example/simple_example.tsv /tmp/legacy_before --plot none --no-plot-tree --events
   ```
2. `git mv medicc2 medicc2-legacy`. **Do not edit its contents** (logger `medicc-main`, version string, options all stay).
3. Update the callers to `medicc2-legacy` (all found with `git grep -n medicc2`):
   * `medicc/test_medicc.py`: the 12 `subprocess.Popen(['python', "medicc2", …])` calls (lines 16–371);
   * `examples/run_all_examples.sh`: the 5 `python medicc2 …` lines;
   * `setup.py`: `scripts=['medicc2-legacy']` for now (Step 46 adds the new `medicc2`);
   * `README.md`: the two usage lines (84, 128). Add one sentence: "The previous OpenFST implementation is available as
     `medicc2-legacy`." Leave the badges and package name `medicc2` alone.
   * `.github/workflows/tests.yml` needs no change: it runs `pytest medicc/test_medicc.py`.

**Test:**
```bash
python medicc2-legacy examples/simple_example/simple_example.tsv /tmp/legacy_after --plot none --no-plot-tree --events
diff -r /tmp/legacy_before /tmp/legacy_after && echo "PASS step 44 (legacy outputs unchanged)"
python -m pytest medicc/test_medicc.py -k help_box -v          # pytest + pytest-timeout installed in medicc-migration
git grep -n -E "\"medicc2\"|python medicc2 " -- medicc examples   # must print nothing
```
(Optional, slow: the full `python -m pytest medicc/test_medicc.py`, like master's CI.)

## Step 45: `medicc/full_fst/cli.py`: the new default `medicc2`, with master's conventions
**REF:** `medicc/full_fst/cli.py`. **Type it with the B2 deviations:**
* parser: keep master names and defaults: `--input-type` default `"t"`, `--input-allele-columns` default
  `"cn_a, cn_b"`, `--plot` default `"auto"`, `-j/--n-cores` default `None` (→ `n_cores or 1`), `--maxcn` without
  `choices` (master's reader raises above 8), and add `-vv/--debug` and `--silent`. New options:
  `--tree-method {nj,fastme}` (default `nj`), `--fastme`, the NNI flags, `--events`, `--chromosomes-bed`, `--regions-bed`.
* logging: **delete `logging.basicConfig`**. Copy master's `medicc2` lines 143–149 (set `logging.getLogger('medicc')`
  to INFO / DEBUG / CRITICAL). Use `logger = logging.getLogger('medicc-main')` with master's `logger.info(...)`
  messages.
* **delete** the "output prefix already exists" refusal (master overwrites).
* before writing trees: copy REF `medicc/core.py::_wrap_tree_for_output` (lines 1010–1036) into cli.py **without its
  `update_branch_lengths` call** (lengths are already final), and **append the MRCA before the normal** (swap the two
  `append` lines), so the order matches master's `(…MRCA…, diploid)`. Apply it to `result.tree`.
* `--tree`: load it with master's `medicc.io.import_tree(args.tree, normal_name=args.normal_name)` (master's rooting
  and checks: unique names, binary, at most one unnamed internal node), then pass the tree object to `reconstruct`.
  User names are kept; see the naming box in Step 36.
* pairwise distances: when `result.pairwise_distances is None` (supplied tree), write a zero DataFrame like master.
* summary: write master's `_summary.tsv` via `medicc.summarize_patient(tree, None, labels, normal_name, events_df)`
  (`pdm` is unused there). If you want REF's json, write it as `<prefix>_full_fst_details.json`.
* end with `logger.info("MEDICC2 has finished.")` instead of `print`.
* plotting: same as REF, but `auto` = bars if fewer than 50 samples, otherwise heatmap (master's rule).

Create the new launcher **`medicc2`** at the repo root (the name is free after Step 44). It is REF's `medicc2`
**without** the env-var and bytecode lines. Its `--version` string is REF's (`… full-fst-l-loh-wgd-2`), so the two
commands are easy to tell apart.
**Test:** `tests/migration/test_45_cli.py`
```python
# <the header>
import os, tempfile
from medicc.full_fst.cli import main
out = tempfile.mkdtemp()
inp = str(ROOT / "examples/simple_example/simple_example.tsv")
main([inp, out, "-j", "2", "--no-plot-tree", "--plot", "none", "--events"])
files = set(os.listdir(out))
for f in ["simple_example_pairwise_distances.tsv", "simple_example_final_tree.new", "simple_example_branch_lengths.tsv",
          "simple_example_final_cn_profiles.tsv", "simple_example_summary.tsv", "simple_example_copynumber_events_df.tsv"]:
    assert f in files, f
from Bio import Phylo
tree = Phylo.read(os.path.join(out, "simple_example_final_tree.new"), "newick")
assert tree.root.name is None and [c.name for c in tree.root.clades][1] == "diploid"   # master shape
clades = lambda t: {n.name: sorted(x.name for x in n.get_terminals()) for n in t.get_nonterminals() if n.name}
assert clades(tree) == {                                          # identical to legacy master medicc2 output
    "internal1": ["taxon_2", "taxon_3"],
    "internal4": ["taxon_2", "taxon_3", "taxon_5"],
    "internal3": ["taxon_1", "taxon_2", "taxon_3", "taxon_5"],
    "internal2": ["taxon_1", "taxon_2", "taxon_3", "taxon_4", "taxon_5"]}
bl = open(os.path.join(out, "simple_example_branch_lengths.tsv")).read().split()
assert {"internal1", "internal2", "internal3", "internal4"} <= set(bl)
main([inp, out, "--no-plot-tree", "--plot", "none"])                          # rerun overwrites (master)
supplied = os.path.join(out, "supplied.new")                                   # --tree keeps user names
# master's import_tree requires every internal node except the root to be named (checked in master)
open(supplied, "w").write("((((taxon_5,(taxon_2,taxon_3)myclade)mid,taxon_1)up,taxon_4)top,diploid);")
out2 = tempfile.mkdtemp()
main([inp, out2, "--tree", supplied, "--no-plot-tree", "--plot", "none"])
t2 = Phylo.read(os.path.join(out2, "simple_example_final_tree.new"), "newick")
assert set(clades(t2)) == {"myclade", "mid", "up", "top"}
print("PASS step 45")
```
Then run both by hand and compare. Both should give tree length 7 on `simple_example` (checked: full model + master NJ = 7, legacy = 7):
```bash
python medicc2        examples/simple_example/simple_example.tsv /tmp/ff -v        # new default (full model)
python medicc2-legacy examples/simple_example/simple_example.tsv /tmp/legacy -v    # master's OpenFST CLI
```

## Step 46: `medicc/full_fst/__init__.py`, `setup.py`, final regression
* `__init__.py`: REF's file **without `FULL_TWO`** (export `Model`, `PairScorer`, `reconstruct`, `distance_matrix`, `Result`).
* `setup.py` (master's file): add `'medicc.full_fst'` to `packages`, set `scripts=['medicc2', 'medicc2-legacy']`, and the native
  assets/libraries to `package_data` (REF's `package_data` globs). Optionally port REF's `BuildEngines` class
  (`ancestors` + `direct_pair` only). Otherwise document `python -B scripts/build.py` as the native build step.
* `doc/medicc2.yml`: add `libboost-headers` and `fastme=2.1.6.3` (bioconda), which the native build and FastME need in CI.
  Add a CI step that runs `python -B scripts/build.py` after `pip install -e .`.
* Port the remaining REF checks if you want them (`tests/check_api.py`, `check_events.py`, `check_capped.py` with
  `reference_pair` removed, and `check_fastme_labels.py`). Adapt `check_api.py` to `internalN` names and to the absence
  of `FULL_TWO`.

**Final test:**
```bash
for t in tests/migration/test_*.py tests/check_pairs.py tests/check_resolve_tie.py tests/check_ancestor_bounds.py; do
  python -B $t || { echo "FAIL $t"; break; }; done
python -m pytest medicc/test_medicc.py -v            # master's suite, now exercising medicc2-legacy
python medicc2 --version && python medicc2-legacy --version   # two distinct commands
```

## Appendix (optional)
* `reference_pair/`: only if you need experiment 20's independent backend. Copy assets and generated headers, type
  `pairwise.cpp` (90 lines; the other headers differ from `direct_pair` by ≤17 lines each), and add it to `build.py`
  and the Step 4 generator list.
* Docs: copy `doc/algorithm_reading_guide.rst`, `doc/full_fst_algorithms.pdf`, and `scripts/build_algorithm_pdf.py` if
  you want them in-repo. Update the guide's `FULL_TWO` / CLI wording to match.
