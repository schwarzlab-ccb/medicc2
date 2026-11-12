# Build plan: a capped full-FST engine for MEDICC2, from master, like a developer

**Updated 2026-11-11.** You build the engine yourself, one file at a time, in the order a developer
with an idea (and no finished code) would build it. Older plans (only on the other machine): `../backups/MIGRATION_PLAN.codex-20261111.md`,
`../backups/MIGRATION_PLAN.claude-v1-20261111.md`.

## <a id="contents"></a>Contents

* [Where this order comes from](#story)
* [How we work](#how)
* [The model in plain words](#plain)
* [Phase 0 — Clean ground](#phase-0) (Steps 0.1–0.4)
* [Phase 1 — An answer key you trust](#phase-1) (Steps 1–3)
* [Phase 2 — Make the model's parts yourself](#phase-2) (Step 4)
* [Phase 3 — The simplest engine that is correct](#phase-3) (Steps 5–11)
* [Phase 4 — Plug it into MEDICC as the default](#phase-4) (Steps 12–21)
* [Phase 5 — Make it fast, without changing one answer](#phase-5) (Steps 22–31)
* [Phase 6 — Tell the story: event histories](#phase-6) (Steps 32–35)
* [Phase 7 — Better trees: NNI and FastME](#phase-7) (Steps 36–40)
* [Phase 8 — Ship it](#phase-8) (Steps 41–46)
* [Appendix: what we keep and what we never build](#appendix)

---

| | Path (relative to `~/Projects/project-medicc`) | Revision | Mamba env |
|---|---|---|---|
| **REF (answer key)** | `medicc-full-fst-lloh-wgd2-capped-x2` | branch `full-fst-lloh-wgd2-capped-x2`, commit `a0ff789` (Codex speed pass on top of tag `openfst-implementation`, `3ece895`); checkout is at `578c60a`, which only adds PGO build files on top | `medicc-codex-ref` |
| **NEW (you build here)** | `feature-full-fst-capped-wgd2` | branch `feature-full-fst-capped-wgd2` (master baseline `04e9bfe`) | `medicc-migration` |

```bash
mamba activate medicc-migration
cd ~/Projects/project-medicc/feature-full-fst-capped-wgd2
```

---

## <a id="story"></a>Where this order comes from

The order follows how the code was **really** built. The records are on `ramses-login` in
`/scratch/cnie/openfst_native_20261111/`:

1. **The engine grew in small steps** (`revisions/development_scripts/`, `DESIGN.rst`, `REPORT.rst`):
   first a 168-line engine (`src/native_unsegmented.cpp`) with plain Dijkstra, no tricks.
   It was checked against an independent oracle. Then, one at a time:
   chromosome splitting → exact score cache → A* → bounds → pair filter.
   After each one: the same answers, and a timing.
   A last speed pass (`performance_20261111/`, REF commit `a0ff789`) removed A* from joins,
   added "prove a known answer", a simpler composition filter, and a bounded cache kept across rows.
2. **Then it went into master** (`master_migration/BRIEF.md`, `steps/`, `patches/`):
   tech stack → native model as **opt-in** (master's default unchanged) → native as **default** → NNI → FastME → FastME fallback.
   We skip the opt-in switch: the engine goes straight in as the default.

Two lessons from those records:
* A speed trick is not always faster. A* alone was **slower** than Dijkstra on most cases (`REPORT.rst`).
  The bounds and the pair filter made it pay off, but only for pair distances. So you measure after every speed step.
* Output differences are of two kinds. A **format** difference (`0.0` vs `0`) is fine.
  A **meaning** difference (a distance, a branch length, a tree) must be explained, never patched over.

---

## <a id="how"></a>How we work

**You are the developer. REF is the answer key, not the source.**

1. The assistant says: **"Please type `<file>` now."** One file per step.
   If a step also needs a one-line hook in another file, the step says so.
2. Each step has: **Problem** (what is wrong or missing now) → **Idea** (how to fix it) → **Type** → **Test** → **Peek**.
3. Try first. Write the test, see it fail, then write the code. Peek at REF only when you are stuck, or after you finish, to compare.
4. Tell the assistant when you are done. The assistant reviews your file. When the test passes, you commit.

**Rules:**

* **Nothing is copied.** You make generated data with your own tools. The only exception is raw input data
  (patient copy numbers such as `meso5.tsv`): you did not create the patients, so you may take that file.
* Use REF's names (`join`, `resolve`, `score`, `ModelConfig`, ...). Then later steps and REF's tests line up.
* Your file does not need to match REF line by line. It must give the same **answers**.
* For a file that exists in master, change only what the step needs.
* **Speed steps never change answers.** Each one re-runs the golden test (Step 11) and the timing (Step 22).

**Lesson tests** live in `tests/migration/`. They are plain scripts, not pytest. Each Python one starts with:

```python
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
```

Run from NEW's root: `python -B tests/migration/test_NN_name.py`.

**C++ lesson tests.** Headers are parts of **one** file (`native.cpp`), in this order:
`types.hpp`, `fst_operations.hpp`, `search.hpp`, `messages.hpp`, `scoring.hpp`
(later also `astar.hpp`, `pair_filter.hpp`, `bounds.hpp`, `history.hpp`). A test includes `<cassert>` and then the list up to its step.

```bash
mkdir -p build/migration
"${CXX:-c++}" -std=c++17 -I. -I"$CONDA_PREFIX/include" tests/migration/test_NN_name.cpp \
  -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lfst -o build/migration/test_NN_name
build/migration/test_NN_name
```

---

## <a id="plain"></a>The model in plain words

* **Profile.** One sample's copy numbers as a string. Each character is one bin, `0`..`8`. `X` separates chromosomes. Example: `"12X3"`.
* **Edge.** One branch of the tree: **parent → child**.
* **A.** Small changes on one edge: losses (`L_LOH`; a loss can go to 0, which is LOH), then gains (`G`). `A = L_LOH ∘ G` (`∘` = "first this, then that").
* **D.** One whole-genome doubling (WGD): every CN doubles, capped at 8: `D(c) = min(2c, 8)`. A WGD costs **1 for the whole genome**.
* **P.** All that can happen on one edge: `P = A + ADA + ADADA` (`+` = "or"). So 0, 1 or 2 WGDs.
* **T(parent, child).** The cheapest cost along one edge. Directed: `T(0→1)` is impossible.
* **K(x, y).** Distance between two samples: `min over z of T(z, x) + T(z, y)`. Example: `K(0, 1) = 1`.
* **Weighted graph (FST).** A compact store for many strings and their costs. OpenFST does the heavy work.
* **Message.** A graph of **all** best ancestor strings for a subtree, plus an `offset` (a cost kept outside the graph).
* **Join** (bottom-up): two children's messages → the parent's message. Keep **every** tied best ancestor.
* **Resolve** (top-down): the parent is fixed → pick one best child string.
* **Bounds.** Lower bound: "it cannot cost less". Upper bound: "we found a way at this cost". They speed up the search. They never change the answer.

**The big idea of the engine:** never build `P` in full (it is huge). For each WGD count `k` (0, 1, 2) build the graph
`(AD)^k ∘ A ∘ child` **lazily** (OpenFST `ComposeFst` makes states only when the search asks), search it, and add `k` once.

---

# <a id="phase-0"></a>Phase 0 — Clean ground

[↑ Contents](#contents)

**Theme:** master runs cleanly on the modern stack, before we add anything.
(The agent's step 1: "never downgrade the new stack; move master forward".)

## Step 0.1: done already

Python 3.12 / NumPy 2 / pandas 3 fixes in `plot.py` (`np.vstack`), `io.py`, `event_reconstruction.py`, packaging pins.
Tests are back in `medicc/`. `pyproject.toml` has `testpaths = ["medicc"]` and skips `tests/migration/`.

```bash
python -m pytest medicc/test_fsts.py medicc/test_medicc.py -q -k 'fstlib or help_box' --timeout=60   # 67 passed
```

## Step 0.2: `medicc/sim.py` — `rcoal` crashes on pandas 3

**Problem.** `python -c "import medicc.sim as s; s.rcoal(4)"` fails with `KeyError: 0`.
`edge.iloc[root][0]` now means "the column named 0", not "the first column".
**Idea.** Use the column names: `['a']`, `['b']`, and `edge.loc[root, 'b']`. Start `reordered` from the first row, not an empty frame.
**Type:** the fix in `rcoal`.
**Test:** `tests/migration/test_00_sim.py`:

```python
import medicc.sim
tree = medicc.sim.rcoal(4)
assert sorted(c.name for c in tree.get_terminals()) == ["diploid", "sp_0001", "sp_0002", "sp_0003", "sp_0004"]
print("PASS step 0.2")
```

**Peek:** `diff -u medicc/sim.py ../medicc-full-fst-lloh-wgd2-capped-x2/medicc/sim.py`.

## Step 0.3: `medicc/plot.py` — old matplotlib calls

**Problem.** The test below fails with `PendingDeprecationWarning: The set_constrained_layout_pads ...`.
Also `groupby('chrom')` on a category column needs `observed=False` in pandas 3.
**Idea.** Use `layout='constrained'` and `fig.get_layout_engine().set(...)`. Add `observed=False` to each `groupby('chrom' ...)`.
**Type:** the fixes in `plot.py`.
**Test:** `tests/migration/test_00_plot.py`:

```python
import warnings
import matplotlib
matplotlib.use("Agg")
import medicc
frame = medicc.io.read_and_parse_input_data(ROOT / "examples/simple_example/simple_example.tsv")
with warnings.catch_warnings():
    warnings.simplefilter("error")
    medicc.plot.plot_cn_profiles(frame)
print("PASS step 0.3")
```

## Step 0.4: `medicc/tools.py` — `iteritems` is gone

**Problem.** `df.iteritems()` does not exist in pandas 2+.
**Idea.** `df.items()`, and `observed=False` in the `groupby`.
**Test:** `grep -n iteritems medicc/tools.py` finds nothing, and `python -c "import medicc.tools"` works.

---

# <a id="phase-1"></a>Phase 1 — An answer key you trust

[↑ Contents](#contents)

**Theme:** before you write a fast engine, have a slow, simple one that you are sure of.
Every later step is checked against it.

## Step 1: `medicc/factory.py` — capped WGD (done)

`create_1step_WGD_fst(..., wgd_x2=True, wgd_cap=True)` caps doubling at 8. Commit `8ec7bc6`.
**Test:** `python -B tests/migration/test_01_wgd_cap.py`.

## Step 2: `tests/oracle.py` — the slow oracle (done)

Builds literal fstlib graphs and pushes small inputs through them. `directed(x, y)` = `T`, `kernel(x, y)` = `K`,
`optimum(left, right)` = all best ancestors.
**Test:** `python -B tests/migration/test_02_oracle.py`.

## Step 3: `tests/make_golden.py` — freeze the oracle's answers

**Problem.** The oracle is slow. You will check the engine many times. You do not want to wait each time.
**Idea.** Run the oracle once on **every** two-bin profile (`"00"` .. `"88"`, 81 words) and save the answers:
`words` (81), `directed` (81×81, `T`) and `undirected` (81×81, `K`). Use `inf` for "impossible".
Save to `medicc/test_data/two_bin_oracle.npz`. (REF has a file with this name and these keys. Yours is made by your own oracle.)
**Type:** `tests/make_golden.py`. Run it once. It can take some minutes.
**Test:** `tests/migration/test_03_golden.py` — it checks the file **makes sense**, not that it matches REF:

```python
import numpy as np
with np.load(ROOT / "medicc/test_data/two_bin_oracle.npz") as f:
    words = f["words"].tolist()
    directed, kernel = f["directed"], f["undirected"]
assert set(words) == {f"{a}{b}" for a in range(9) for b in range(9)}
assert directed.shape == kernel.shape == (81, 81)
assert np.array_equal(kernel, kernel.T)                       # K is symmetric
assert np.all(np.diag(directed) == 0)                         # no change costs 0
assert np.isinf(directed[words.index("00"), words.index("11")])  # cannot gain from 0
assert np.all(kernel <= directed)                            # the child itself can be the ancestor
print("PASS step 3")
```

**Extra check (one time):** compare your file with REF's. They should be equal. If not, one of the two oracles is wrong. Find out which.

```bash
python -c "import numpy as np; a=np.load('medicc/test_data/two_bin_oracle.npz'); b=np.load('../medicc-full-fst-lloh-wgd2-capped-x2/medicc/test_data/two_bin_oracle.npz'); print(all(np.array_equal(a[k], b[k]) for k in ('directed','undirected')))"
```

Delete the old `tests/migration/test_03_native_assets.py` (it asked you to copy files).

---

# <a id="phase-2"></a>Phase 2 — Make the model's parts yourself

[↑ Contents](#contents)

**Theme:** the C++ engine needs the small transducers `L_LOH`, `G`, `A` and `AD`. You already know how to build them in Python. Build them, and write them to a file that C++ can read.

## Step 4: `tools/export_model_parts.py` — write `medicc/objects/ancestral_parts.txt`

**Problem.** C++ cannot call fstlib. The engine needs the primitives as plain data.
**Idea.** Build the parts with `medicc.factory` (as `tests/oracle.py::get_parts` does):

* `L_LOH` and `G` from the factory.
* `A = L_LOH ∘ G`, then make it small (determinize/minimize where fstlib allows).
* `D` = capped doubling **with cost 0** (the WGD cost is added once per genome by the engine, not per bin). See `doubling()` in your oracle.
* `AD = A ∘ D`, made small.

Write each one as a text record:

```
NAME n_states start n_arcs
final cost of state 0 ... state n-1          (one line)
source ilabel olabel cost target             (one line per arc)
```

Labels in the file: CN `0..8` → `0..8`, `X` → `9`. (C++ adds 1 later, so OpenFST label 0 can stay "epsilon".)
Choose and document a number for "not final".
**Type:** `tools/export_model_parts.py`. Run it.
**Test:** `tests/migration/test_04_model_parts.py` — reads your file back:

```python
tokens = iter((ROOT / "medicc/objects/ancestral_parts.txt").read_text().split())
seen = set()
for name in tokens:
    seen.add(name)
    states, start, count = (int(next(tokens)) for _ in range(3))
    assert states > 0 and 0 <= start < states
    finals = [float(next(tokens)) for _ in range(states)]
    for _ in range(count):
        s, i, o, w, t = (int(next(tokens)) for _ in range(5))
        assert 0 <= s < states and 0 <= t < states and 0 <= i <= 9 and 0 <= o <= 9 and w >= 0
assert {"L_LOH", "G", "A", "AD"} <= seen
print("PASS step 4")
```

Your state and arc counts may differ from REF's file (REF: `A` 212 states / 15688 arcs, `AD` 74 / 3108). That is fine.
The real proof is Step 11: the engine with **your** file must match the golden answers.
**Peek:** how REF's loader reads the file: `load_model` in `REF/medicc/native/src/native.cpp`.

---

# <a id="phase-3"></a>Phase 3 — The simplest engine that is correct

[↑ Contents](#contents)

**Theme:** a small C++ engine with plain Dijkstra search. No splitting, no cache, no A*, no bounds.
It is slow but correct. (This is the real first version: `src/native_unsegmented.cpp`, 168 lines.)
All files are in `medicc/native/src/`. Read [the model in plain words](#plain) first.

## Step 5: `types.hpp` — the shared words

**Problem.** Every C++ file needs the same names for a cost, a graph, a message and the model.
**Idea.** `W` = `fst::TropicalWeightTpl<double>` (costs add along a path; take the minimum between paths).
`Arc`, `Fst`, `Vec`, `Ptr` = OpenFST graph types. `Message {Vec fst; double offset; int length;}`.
`Model {Vec a, ad; int max_wgds = 2;}` — **only what you need now**. Later steps add fields.
**Type:** `types.hpp`.
**Test:** rename the existing test: `git mv tests/migration/test_04_types.cpp tests/migration/test_05_types.cpp`. Compile and run it.

## Step 6: `fst_operations.hpp` — strings become graphs

**Problem.** The engine gets strings like `"1X2"`. OpenFST needs graphs and integer labels.
**Idea.** `label('0') = 1 … label('8') = 9`, `label('X') = 10` (0 stays "epsilon"). `symbol` goes back.
`linear("1X2")` = a graph that accepts only that string, at cost 0. `compatible(x, y)` throws if lengths or `X` positions differ.
`compose(a, b)` = a **lazy** `ComposeFst` (sort the right side's arcs with `ArcSortFst` first). `check(f)` throws on an OpenFST error.
**Type:** `fst_operations.hpp`.
**Test:** `tests/migration/test_06_fst_operations.cpp`:

```cpp
int main() {
    assert(label('0') == 1 && label('X') == 10 && symbol(9) == '8');
    Vec f = linear("1X2");
    assert(f.NumStates() == 4);
    assert(fst::ShortestDistance(f).Value() == 0);
    bool rejected = false;
    try { compatible("1X2", "12X"); } catch (const std::runtime_error &) { rejected = true; }
    assert(rejected);
}
```

## Step 7: `search.hpp` — shortest distance that keeps all ties

**Problem.** You need the best cost through a lazy graph. For ancestors you need **every** best path, not one.
And you want to stop early when a search cannot beat a known cost (`cap`).
**Idea.** Use OpenFST's `ShortestDistance` with your own queue (`BoundQueue`, based on a shortest-first queue).
The queue says "empty" when the next state costs **more than** `cap`. (Strictly more: equal states hold more tied answers.)
After the search, keep each arc with `d(u) + w == d(v)` and each final with `d(u) + final == best`. `Connect` removes dead ends.
That graph holds **all** best paths. Return `Search {double cost; Vec tight;}`.
Also: `canonical(f)` (clean up: remove epsilons, determinize, minimize) and `choose(f)` (one shortest string).
**Type:** `search.hpp`.
**Test:** `tests/migration/test_07_search.cpp`:

```cpp
int main() {
    Vec f = linear("1"), other = linear("2");
    fst::Union(&f, other);
    auto answer = search(f, infinity, true);
    assert(answer.cost == 0);
    auto accepted = [&](const std::string &word) {
        return std::isfinite(fst::ShortestDistance(*compose(linear(word), answer.tight)).Value());
    };
    assert(accepted("1") && accepted("2") && !accepted("3"));   // both ties kept
    assert(!std::isfinite(search(f, -1, true).cost));           // cap too small
}
```

**Peek:** `REF/medicc/native/src/search.hpp`. REF's queue also has an A* part (`Priority`, `Heuristic`). Leave that out for now; it comes in Step 27.

## Step 8: `messages.hpp` — one edge, join and resolve

**Problem.** For a tree you need two operations: bottom-up **join** and top-down **resolve**.
**Idea.**
* `pull(model, child, k)` = the graph of all parents that can reach `child` with `k` WGDs:
  `A ∘ child`, then `AD ∘ …` `k` times, then project to the input side. All lazy.
* `join(model, left, right)`: for each `(k, j)` in `0..2 × 0..2`, compose `pull(left, k)` with `pull(right, j)`,
  search, add `k + j`. Keep the minimum. **Union** the tight graphs of every `(k, j)` that ties. Store cost in `offset`.
* `resolve(model, parent, child)`: for each `k`, compose `linear(parent)` with the edge graph, find the best child string; add `k` and the child's offset.

**Type:** `messages.hpp`.
**Test:** waits for Step 11 (it needs the model file loaded).
**Peek:** REF's `join`/`resolve` also call `split` and bounds. Leave those out. In REF the graph builder is `make_graph` in `astar.hpp`. You move yours there in Step 27.

## Step 9: `scoring.hpp` — the pair distance

**Problem.** The tree builder needs `T(x, y)` and `K(x, y)` for profile strings.
**Idea.** `fixed(model, x, y, k, j, directed)` = one search for one WGD choice. `score(model, x, y, directed)` = minimum over all choices, plus WGD costs.
Directed: `x` is the parent, so only the child side has edges. Symmetric: both sides, as in `join`.
**Type:** `scoring.hpp`.
**Test:** waits for Step 11.

## Step 10: `native.cpp` — load the model, open a C door

**Problem.** Python cannot call C++ functions directly. And the model file must be read.
**Idea.** `#include` all headers once. `load_model(path)` reads your Step 4 file (add 1 to each label).
Plain C functions (`extern "C"`): `nf_load`, `nf_score`, `nf_join`, `nf_resolve`, `am_string`, `am_drop`, `am_free`,
`am_offset`, `am_error`. Never throw across the C door: catch, store the text in a `thread_local` string, return null or −1.
**Type:** `native.cpp`. Build it by hand for now:

```bash
"${CXX:-c++}" -std=c++17 -O2 -shared -fPIC -I"$CONDA_PREFIX/include" medicc/native/src/native.cpp \
  -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lfst -o medicc/native/_engine.so
```

**Test:** `tests/migration/test_10_native.cpp` (include `native.cpp` itself):

```cpp
#include <cassert>
#include "medicc/native/src/native.cpp"
int main() {
    void *model = nf_load("medicc/objects/ancestral_parts.txt", 0);
    assert(model);
    assert(nf_score(model, "1X1", "2X2", 1) == 1);          // one WGD for the whole genome
    assert(!std::isfinite(nf_score(model, "0", "1", 1)));
    assert(nf_score(model, "0", "1", 0) == 1);
    void *a = am_string("1"), *b = am_string("2");
    void *joined = nf_join(model, a, b);
    assert(joined && am_offset(joined) == 1);
    am_drop(joined); am_drop(a); am_drop(b); am_free(model);
}
```

## Step 11: `medicc/native/model.py` — Python wrapper, and the big check

**Problem.** You want to call the engine from Python, and compare it with the oracle on all 81×81 pairs.
**Idea.** `ctypes`: find `_engine*.so` next to the file, declare each function's argument and return types.
`ModelConfig(no_wgd=False, total_copy_numbers=False)` (a small frozen dataclass).
`Model`: owns the native pointer; `score`, `string`, `join`, `resolve`; `close()` and `with` support.
`Message`: owns one native message; `close()`. Turn a null/−1 return into a Python exception with the error text.
**Type:** `touch medicc/native/__init__.py`, then `medicc/native/model.py`.
**Test:** `tests/migration/test_11_engine_vs_oracle.py` — **the milestone of this phase**:

```python
import numpy as np
from medicc.native.model import Model, ModelConfig
with np.load(ROOT / "medicc/test_data/two_bin_oracle.npz") as f:
    words, directed, kernel = f["words"].tolist(), f["directed"], f["undirected"]
with Model() as model:
    for i, x in enumerate(words):
        for j, y in enumerate(words):
            assert model.score(x, y, directed=True) == directed[i, j], (x, y)
            assert model.score(x, y) == kernel[i, j], (x, y)
    a, b = model.string("11"), model.string("22")
    joined = model.join(a, b)
    assert model.resolve("11", joined)[1] == model.resolve("22", joined)[1]  # both ties kept
    for m in (joined, a, b):
        m.close()
with Model(ModelConfig(no_wgd=True)) as model:
    assert model.score("1X1", "2X2", directed=True) == 2
print("PASS step 11: engine == oracle on all two-bin pairs")
```

If it fails: is it your model file (Step 4) or your engine? Try REF's `ancestral_parts.txt` **once, only to find out which**. Then fix your own.
Then write down how long the test takes. That is your first timing.

---

# <a id="phase-4"></a>Phase 4 — Plug it into MEDICC as the default

[↑ Contents](#contents)

**Theme:** MEDICC uses your engine by default. A custom `--fst PATH` still uses master's fstlib code.
There is no `--native` switch: we change master directly.

The engine is still slow here (Phase 5 makes it fast). So in this phase, run the small examples and the fast test groups.
Run the full suite after Phase 5.

## Step 12: `setup.py` — build the engine with the package

**Problem.** The hand build command (Step 10) is easy to forget.
**Idea.** Add a third `Extension("medicc.native._engine", ["medicc/native/src/native.cpp"], ...)` next to the two fstlib ones.
C++17, `-O3`, link `fst`, same include/lib/rpath handling. It is loaded by ctypes, not imported.
**Test:** `rm -f medicc/native/_engine*.so && python setup.py build_ext --inplace && python -B tests/migration/test_11_engine_vs_oracle.py`.

## Step 13: `pyproject.toml` — ship the new files

**Problem.** `pip install .` would leave out `medicc.native`, the `.hpp` files and `ancestral_parts.txt`.
**Idea.** Add `medicc.native` to packages; package data `objects/*.txt` and `"medicc.native" = ["src/*.hpp", "src/*.cpp"]`. Keep the pytest block.
**Test:** `python -c "import tomllib; d=tomllib.load(open('pyproject.toml','rb')); assert 'medicc.native' in d['tool']['setuptools']['packages']"`.

## Step 14: `medicc/native/scheduler.py` — many pairs, many threads

**Problem.** A tree needs all pair distances, and one join per inner node. One thread is slow.
**Idea.** `pairwise_distances(cn_str_dict, config, n_cores)` → a symmetric `DataFrame` in sample order.
`reconstruct(children, leaves, root, normal, config, n_cores)`: bottom-up joins when both children are ready, then top-down resolves.
**Each thread has its own `Model`** (a model has caches; it is not thread-safe). Close all models at the end.
**Type:** `scheduler.py`. Then put the exports in `medicc/native/__init__.py`.
**Test:** `tests/migration/test_14_scheduler.py`:

```python
import pandas as pd
from medicc.native import ModelConfig, pairwise_distances, reconstruct
samples = {"b": "22", "a": "11"}
one = pairwise_distances(samples, ModelConfig(), n_cores=1)
two = pairwise_distances(samples, ModelConfig(), n_cores=2)
pd.testing.assert_frame_equal(one, two)
assert list(one.index) == ["b", "a"] and one.loc["b", "a"] == 1
results = [reconstruct({"r": ["a", "b"]}, samples, "r", "11", ModelConfig(), n_cores=n) for n in (1, 2)]
assert results[0] == results[1] == {"b": "22", "a": "11", "r": "11"}
print("PASS step 14")
```

## Step 15: `medicc/io.py` — the default model is native

**Problem.** `read_fst()` builds master's old fstlib model.
**Idea.** `read_fst(no_wgd=False, total_copy_numbers=False, user_fst=None)`: with no `user_fst` it returns `ModelConfig(no_wgd, total_copy_numbers)`.
With `user_fst=PATH` it loads that fstlib FST, as master does. Allow a single tumour sample for the native model (`validate_input(..., min_samples=...)`).
`write_branch_lengths` leaves out the normal.
**Test:** `tests/migration/test_15_io.py`:

```python
import fstlib
from medicc import io
from medicc.native import ModelConfig
assert io.read_fst() == ModelConfig()
assert io.read_fst(no_wgd=True, total_copy_numbers=True) == ModelConfig(True, True)
assert isinstance(io.read_fst(user_fst=str(ROOT / "medicc/objects/wgd_asymm.fst")), fstlib.Fst)
print("PASS step 15")
```

## Step 16: `medicc/ancestors.py` — send native work to the scheduler

**Idea.** If `fst` is a `ModelConfig`: turn the tree into `children/leaves/root/normal`, call `reconstruct`. Otherwise master's FSA code. Add `n_cores`.
**Test:** `tests/migration/test_16_ancestors.py`:

```python
from io import StringIO
from Bio import Phylo
from medicc.ancestors import reconstruct_ancestors
from medicc.native import ModelConfig
tree = Phylo.read(StringIO("(diploid,a,b)r;"), "newick")
result = reconstruct_ancestors(tree, {"diploid": "11", "a": "11", "b": "22"}, ModelConfig(), "diploid", n_cores=2)
assert result["r"] == "11" and result["b"] == "22"
print("PASS step 16")
```

## Step 17: `medicc/core.py` — the pipeline uses strings for the native model

**Idea.** In `main`: if the model is a `ModelConfig`, use strings (`create_cn_str_dict`), native distances, native reconstruction, native branch lengths (`reconstructed_tree_length`), and turn strings back into a table (`create_df_from_cn_str_dict`). Reject `prune_weight > 0` (`_check_native_options`). Keep the six return values. Keep master's FSA path for a custom FST.
Type it in two turns: (a) the string/table helpers, (b) `main`.
**Test:** `tests/migration/test_17_core.py` — take `meso5.tsv` and `meso5.new` from REF (raw input data, allowed):

```python
import medicc
frame = medicc.io.read_and_parse_input_data(ROOT / "medicc/test_data/meso5.tsv")
tree = medicc.io.import_tree(ROOT / "medicc/test_data/meso5.new")
out = medicc.main(frame, medicc.io.read_fst(), input_tree=tree, n_cores=2, reconstruct_events=False)
assert out[3].total_branch_length() == 67
print("PASS step 17")
```

## Step 18: `medicc/__init__.py` — import `medicc.native`

One line. **Test:** `python -c "import medicc; medicc.native.ModelConfig()"`.

## Step 19: `medicc2` — the CLI uses the native model

**Idea.** By default call `read_fst(no_wgd=..., total_copy_numbers=...)`. `--fst PATH` = custom fstlib model.
`--prune-weight > 0` with the native model fails early with a clear message. Keep all other options.
**Test:**

```bash
python ./medicc2 examples/simple_example/simple_example.tsv /tmp/native_out --plot none --no-plot-tree
python ./medicc2 examples/simple_example/simple_example.tsv /tmp/fst_out --fst medicc/objects/wgd_asymm.fst --plot none --no-plot-tree
```

Both run. The second one must give master's old numbers.

## Step 20: `medicc/test_fsts.py` — pin the old numbers to the old model

**Problem.** These tests check master's old model numbers. They now get the native model by default.
**Idea.** Load the old model **explicitly** (`read_fst(user_fst=".../wgd_asymm.fst")`). **Never** change an expected number here.
**Test:** `python -m pytest medicc/test_fsts.py -q`.

## Step 21: `medicc/test_medicc.py` — expected changes, listed

**Problem.** Some integration tests check old default outputs.
**Idea.** For each test that fails: is it a **format** difference (accept) or a **meaning** difference (the new model, on purpose)?
Compare the new number with REF's result for the same input (same model, so it should match).
Write each changed number in the commit message: test name, old value, new value, why.
**Test:** `python -m pytest medicc/test_medicc.py -q -k 'not gundem and not ipynb' --timeout=900`.

**Milestone:** MEDICC runs your engine by default, on real data.

---

# <a id="phase-5"></a>Phase 5 — Make it fast, without changing one answer

[↑ Contents](#contents)

**Theme:** add speed tricks in the order they were really added. After **each** step:
1. `python -B tests/migration/test_11_engine_vs_oracle.py` passes (same answers).
2. `python -B tests/migration/bench.py` — write down the time. Did it get faster? If not, say why.

**Where A\* is used (final design, REF commit `a0ff789`):**
* **Pair distances** (`score`, not directed): A\* with suffix potentials, inside the pair filter (Step 27).
* **Joins, resolves, directed scores:** plain bounded Dijkstra. No A\*.

The real history first tried A\* for **all** searches, with a guess from a small "quotient" graph.
Building that small graph cost more time than it saved. Codex removed it in `a0ff789`. You skip it.

## Step 22: `tests/migration/bench.py` — measure first

**Problem.** You cannot say "faster" without a number.
**Idea.** Time three things on a real example with 1 core: the full pairwise distance matrix, one fixed-tree reconstruction (`meso5`), and the slowest single pair you find. Print the times.
Save the results of the first run (distances, tree length) and check that later runs give the **same** numbers.
**Test:** it runs and prints times. Commit it with your first numbers in the commit message.

## Step 23: `messages.hpp` (+ `fst_operations.hpp`, `scoring.hpp`) — split at chromosomes

**Problem.** A long profile gives one huge search. But `A` and `AD` restart at each `X`.
**Idea.** For a **fixed** WGD choice, costs add over chromosomes. So search each chromosome alone and add.
Never choose the WGD count per chromosome: one WGD count for the whole genome.
For messages: cut after an `X` only when **both** messages have exactly one state there (`layers`, `slice`, `split`). Otherwise the pieces are linked.
Score the shortest chromosomes first: a cheap piece fills the budget early and stops hopeless WGD choices sooner.
**Hooks:** `layers` goes in `fst_operations.hpp`; `scoring.hpp::score` splits strings at `X`.
**Test:** `test_11` + bench, plus `tests/migration/test_23_split.cpp`:

```cpp
int main() {
    Message a{linear("1X2"), 0, 3}, b{linear("2X3"), 0, 3};
    auto blocks = split(a, b);
    assert(blocks.size() == 2);
    assert(blocks[0].first.length == 2 && blocks[1].first.length == 1);
}
```

Also test empty chromosomes: `"1XX2"`.

## Step 24: `bounds.hpp` (new, small) — remember exact pair scores

**Problem.** The same chromosome pairs come back again and again in a distance matrix.
**Idea.** A cache in `Model`: key = `(x, y, k, j, directed)` → cost (`score_key`, `cached_score`). `nf_clear` empties it.
**Careful:** store only a **finished** result. A search that stopped at `cap` did not prove "impossible".
For now the cache has no size limit, and the scheduler clears it after each row. Step 31 changes this.
**Test:** `test_11` + bench. Add a test that scores with a small `cap` first, then without, and gets the right answer.

## Step 25: `tools/generate_chromosome_bound.py` — a table of boundary costs

**Problem.** For bounds you need the cheapest possible cost between two neighbouring columns, for each WGD pair `(k_left, k_right)` (9 pairs).
**Idea.** Compute it once in Python (min-plus over the copy-number transitions). Write it as C++ constants into `medicc/native/src/chromosome_constants.hpp`. The WGD charge is **not** in the table.
**Test:** run it twice; the output is the same. Then compare once with REF's header: `cmp medicc/native/src/chromosome_constants.hpp ../medicc-full-fst-lloh-wgd2-capped-x2/medicc/native/src/chromosome_constants.hpp`.

## Step 26: `bounds.hpp` — upper and lower bounds

**Problem.** The search does not know early that a WGD choice is hopeless.
**Idea.**
* Upper bound: build **one** valid history by hand (`changes`, `doubled`, `witness`). Its cost is a cap for the search.
* Lower bound: the boundary table (`string_lower`, `language_lower`).
* `branch_order`: try WGD choices from the lowest lower bound. Skip a choice whose lower bound is worse than the best found. In a **join**, skip only **strictly** worse choices (keep ties).
* If lower bound = upper bound, the answer is known. Cache it and do not search.

**Test:** `test_11` + bench, plus `tests/migration/test_26_bounds.cpp`:

```cpp
int main() {
    assert(doubled("5X3", 1) == "8X6");
    assert(changes("22", "11") == 1);
    assert(changes("2X2", "1X1") == 2);
    assert(!std::isfinite(changes("0", "1")));
}
```

## Step 27: `pair_filter.hpp` + `astar.hpp` — A\* for pairs, and drop arcs that are too expensive

**Problem.** For a pair, OpenFST still builds composition states that cannot fit the cap. Dijkstra also looks at states that cannot lead to a cheap answer.
**Idea.**
* **Suffix potentials.** For each CN string, compose it with `A` once (or with `a_even`: after a WGD every CN is even). Compute forward and backward distances on that small graph. The backward distance is a lower bound on "cost left" (`SuffixData`, cached per string in `suffix_cache`).
* **A\*.** The search orders states by `g + h`: cost so far plus the guess of cost left. `h` = left suffix + right suffix. The guess is never too high, so the answer stays exact.
* **Filter.** A custom `ComposeFilter` (`CostComposeFilter`) sees the current cap (`ComposeBudget`). It drops an arc when prefix + arc + suffix is above the cap, **before** OpenFST makes the state. `make_pair_graph` builds the filtered graph.
* Only for pair scores. Joins, resolves and directed scores stay plain Dijkstra.

**Hooks:**
* `search.hpp`: the queue takes an optional `Heuristic*` (`Priority` = `g + h`) and an optional `ComposeBudget*`. Null = Dijkstra, as before.
* `astar.hpp` (new, small): `Heuristic` (a function `evaluate` plus a cache of its values). Move your graph builder from `messages.hpp` here as `make_graph` (only `compose`, no heuristic).
* `Model` gets `a_even` and `suffix_cache`. `scoring.hpp::fixed` uses `make_pair_graph` when not directed.

**Test:** `test_11` + bench, plus `tests/migration/test_27_astar.cpp`:

```cpp
int main() {
    Heuristic h;
    int calls = 0;
    h.evaluate = [&](int s) { ++calls; return double(s); };
    assert(h(2) == 2 && h(2) == 2 && calls == 1);   // guesses are cached
}
```

This is where the big speed-up came in the real history.

## Step 28: `medicc/native/shortening.py` — drop columns that never change

**Problem.** In a pair, many columns are the same as their neighbour in **both** profiles.
**Idea.** `prepare(profile)` marks where each profile changes (once per sample). `pair(left, right)` keeps only columns where **either** changes, and always keeps `X`. Only for pair distances. Hook: the scheduler uses it.
**Test:** `tests/migration/test_28_shortening.py`:

```python
from medicc.native.shortening import prepare, pair
from medicc.native.model import Model
x, y = "111X22", "122X44"
assert pair(prepare(x), prepare(y)) == ("11X2", "12X4")
assert pair(prepare("1XX1"), prepare("2XX2")) == ("1XX1", "2XX2")
with Model() as model:
    assert model.score(x, y) == model.score(*pair(prepare(x), prepare(y)))
print("PASS step 28")
```

## Step 29: `scoring.hpp` — prove a known answer instead of finding it again

**Problem.** Often the upper bound from Step 26 is already the true answer. The search still works hard to find that same cost again.
**Idea.** Costs are whole numbers. So ask a smaller question: "Is there a path that is **strictly cheaper**?"
* In `score`: search with `room = min(best - other - 1, upper - 1)`.
  If the search finds nothing and `room == upper - 1`, then `upper` is the exact answer. Cache it.
* In `fixed` (pairs only): `greedy_pair` walks the filtered graph once, always taking the arc with the smallest `weight + h`.
  It gives one complete path = a feasible cost. Set `cap = min(cap, feasible - 1)` and search only for cheaper paths. Answer = `min(feasible, search result)`.
* **Careful:** the greedy cost is an upper bound, never the answer by itself. Do not set `active_g` from a greedy path: a greedy prefix is not a shortest path.

**Test:** `test_11` + bench, plus `tests/migration/test_29_certify.cpp`: load the model, check `score("22", "11", false) == 1`, and check that `fixed(..., cap = 0)` returns infinity and is **not** cached.

## Step 30: `pair_filter.hpp` — a simpler composition filter

**Problem.** OpenFST's default composition filter tracks epsilon (empty-label) moves. Pair graphs have **no** epsilon arcs: every arc reads one symbol on each side. That tracking is wasted work.
**Idea.** Base `CostComposeFilter` on `fst::NullComposeFilter<fst::SortedMatcher<Fst>>`, with a matching `GenericComposeStateTable` (`PairMatcher`, `PairFilter`, `PairStateTable`). Only the pair graph changes.
**Test:** `test_11` + bench. Small change, small gain. Write down the number.

## Step 31: `bounds.hpp` + `scheduler.py` — keep the score cache across rows, with a size limit

**Problem.** Step 24 clears the cache after each matrix row. But the next row asks many of the same chromosome pairs.
If you keep it forever, memory grows without limit.
**Idea.**
* `cache_score(m, key, value)`: at most 65536 entries and 8 MiB of key text (`score_key_bytes`). When full, clear all and start again. Clearing only forgets exact answers; it never changes a result.
* `nf_finish_row`: clear only `suffix_cache` (the big per-string FSTs). Keep `chromosome_scores`. `nf_clear` still clears everything.
* `Model._finish_row()` in `model.py`. The scheduler calls it after each row instead of `clear()`.

**Test:** `test_11` + bench, plus `tests/migration/test_31_cache.py`:

```python
from medicc.native.model import Model
with Model() as model:
    first = model.score("1221X3", "2112X1")
    model._finish_row()
    assert model.score("1221X3", "2112X1") == first
print("PASS step 31")
```

**Phase check:** compare your bench numbers from Step 22 and now.
Then run the full suite once: `python -m pytest medicc/test_fsts.py medicc/test_medicc.py -q -k 'not ipynb' --timeout=900` (this includes the slow `gundem` tests).

---

# <a id="phase-6"></a>Phase 6 — Tell the story: event histories

[↑ Contents](#contents)

**Theme:** a distance is a number. Users also want to know **what happened**: which losses, gains and WGDs, where.

## Step 32: `history.hpp` — one concrete history for an edge

**Problem.** The search knows the best cost, but not the steps.
**Idea.** Pick a best WGD count. Walk the stages (loss, gain, doubling) and record the profile after each one, using OpenFST composition state tuples (`trace_stages`). Output: the score, then lines `type<TAB>before<TAB>after`. Keep empty chromosomes (`chromosome_words("1XX2X") == {"1","","2",""}`).
`Model` gets `loss`, `gain`, `doubling`. Add `nf_history` in `native.cpp`.
**Test:** `tests/migration/test_32_history.cpp` with `chromosome_words`, and a call to `nf_history`.

## Step 33: `medicc/native/model.py` — `Model.history`

Returns `(score, stages)`, each stage a dict `type, before, after`.
**Test:**

```python
from medicc.native.model import Model
with Model() as model:
    cost, stages = model.history("1X1", "2X2")
    current = "1X1"
    for stage in stages:
        assert stage["before"] == current
        current = stage["after"]
    assert current == "2X2" and cost == 1
    assert sum(s["type"] == "wgd" for s in stages) == 1
print("PASS step 33")
```

## Step 34: `medicc/event_reconstruction.py` — stages become events

**Idea.** Split each loss/gain stage into unit interval events (`_unit_intervals`, `_native_stage_events`). A WGD stage = one genome-wide event.
Check: replaying the events gives the child, and the number of events = branch length.
Also fix: two WGDs with the same index must both stay in the BED overlap output.
**Test:**

```python
import medicc
frame = medicc.io.read_and_parse_input_data(ROOT / "medicc/test_data/meso5.tsv")
tree = medicc.io.import_tree(ROOT / "medicc/test_data/meso5.new")
_, _, _, tree, _, events = medicc.main(frame, medicc.io.read_fst(), input_tree=tree, reconstruct_events=True)
assert len(events) == tree.total_branch_length() == 67
print("PASS step 34")
```

## Step 35: `medicc/bootstrap.py` — same settings in every resample

Each run uses the same model, no-WGD choice, tree method and normal name. Fewer than 5 runs is an error.
**Test:** `run_bootstrap(frame, original_tree=tree, N_bootstrap=5, seed=1, fst=config, n_cores=2, show_progress=False)` returns 5 counted runs.

---

# <a id="phase-7"></a>Phase 7 — Better trees: NNI and FastME

[↑ Contents](#contents)

**Theme:** the engine gives good distances. Now get better tree shapes from them.

## Step 36: `medicc/nni.py` — try nearby trees

**Problem.** NJ gives one tree. A neighbour tree can be shorter.
**Idea.** `nni_neighbours(tree)`: at each inner edge, swap a subtree across it (2 moves per edge); the normal stays the outgroup; return copies.
`topology_key(tree)`: same key for the same shape, whatever the inner names or child order.
`nni_search(tree, evaluate, max_iter)`: steepest descent; take the best **strict** improvement; skip seen shapes; return `(best_tree, trace)`.
It knows nothing about the engine: you pass `evaluate`.
**Test:**

```python
from io import StringIO
from Bio import Phylo
from medicc.nni import nni_neighbours, topology_key, nni_search
tree = Phylo.read(StringIO("(diploid,(a,b)i1,(c,d)i2)i0;"), "newick")
original = topology_key(tree)
neighbours = list(nni_neighbours(tree))
assert len(neighbours) == 4 and len({topology_key(t) for t in neighbours}) == 4
assert topology_key(tree) == original
best, trace = nni_search(tree, lambda _: 7)
assert topology_key(best) == original and len(trace) == 1
print("PASS step 36")
```

## Step 37: `medicc/core.py` + `medicc2` — wire NNI in

`main(..., nni_mode=..., nni_max_iter=...)`; `evaluate` = reconstructed tree length. CLI `--nni-mode`, `--nni-max-iter`. Write `PREFIX_nni_trace.tsv`.
**Test:** on `meso5` with a bad start tree the length goes down (REF: 98 → 67). Same result with `-j 1` and `-j 2`.

## Step 38: `medicc/core.py` — FastME

**Idea.** `infer_tree_topology(..., method='fastme')`: write a PHYLIP matrix to a temp folder, run `fastme -i D.phy -o T.nwk -m B -n -s -T 1 -z 104000`, read the tree, root on the normal like NJ does.
Long sample names break FastME (64-byte buffer): use safe temporary names and map them back.
**Test:** a sample name of 100 characters comes back unchanged.

## Step 39: `medicc/core.py` — fall back to NJ

`select_tree_method`: use NJ when asked (`--tree-method nj`), when FastME is not found, or with fewer than 4 samples. Log why.
**Test:** with `fastme="no-such-program"` the run finishes with NJ and logs a warning.

## Step 40: `medicc/conftest.py` + `medicc2` — FastME in tests and CLI

`conftest.py` finds FastME (or `FASTME=...`) and skips only tests marked as needing it. CLI: `--tree-method {fastme,nj}`, `--fastme PATH`.
**Test:** `python -m pytest medicc --collect-only -q`.

---

# <a id="phase-8"></a>Phase 8 — Ship it

[↑ Contents](#contents)

**Theme:** lesson tests explain pieces. Now real tests check the whole thing, and others can install it.

| Step | File | What it checks |
|---|---|---|
| 41 | `medicc/test_native.py` | Your golden fixture, correlated joins, errors after `close()`, no-WGD, event replay, CLI. Turn your best lesson tests into pytest tests. |
| 42 | `medicc/test_nni.py` | Neighbour shapes; tree length goes down on a fixed example |
| 43 | `medicc/test_fastme.py` | Name round trips, missing FastME, explicit NJ, small samples |
| 44 | `.github/workflows/tests.yml` | CI builds the package and runs the old groups **plus** native, NNI and FastME tests |
| 45 | `doc/medicc2.yml` | Add FastME. Keep Python 3.12 / OpenFST 1.8.4 |
| 46 | `README.md`, `CHANGELOG.md` | Native vs `--fst`, no positive pruning, total CN, NNI, FastME, build command. One CHANGELOG entry per phase |

Write your tests first. **Then** read REF's test files and add any case you missed.

**Final check:**

```bash
rm -f medicc/native/_engine*.so && python setup.py build_ext --inplace
MPLBACKEND=Agg python -m pytest medicc -k 'not ipynb' --timeout=600 --durations=15
for t in tests/migration/test_*.py; do python -B "$t" || break; done
```

Then run the same examples in REF and NEW, with a fixed tree and without. Compare distances, branch lengths, tree length and event counts **as numbers**.
Ties may pick a different ancestor string; that is fine if the cost is the same and the events replay.
Last: install the package and run `import medicc` and `medicc2 --help` from **outside** the checkout (this finds missing data files).

---

## <a id="appendix"></a>Appendix: what we keep and what we never build

**Keep:** `fstlib` submodule (rev `8216876`) and symlink, `nj.py`, `stats.py`, `factory.py` (with Step 1), master's `objects/*.fst` and `*.bed`,
`tools/medicc-*.py`, examples, notebooks, `tests/oracle.py`, `tests/migration/`.

**Never build** (from older plans or experiments; not needed):
`medicc/full_fst`, separate `direct_pair/`/`ancestors/` engines, Boost.Graph, `scripts/build.py`, `FULL_TWO`, upper/lower FST pairs,
`tree_hash.py`, the quotient A* heuristic (`quotient`, `child_lattice`, `Pattern`; removed in `a0ff789`), SPR, simulated annealing, plateau/frontier NNI, `--search`/`--diagnostics`/`NF_*` switches, stats counters, `medicc2-legacy`.
REF's `medicc/cext/full_fst/.../build/*.so` files are old build leftovers.
