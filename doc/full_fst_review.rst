Reviewing the full L_LOH / two-WGD reference
==========================================

Capped-WGD branch
-----------------

This experimental branch changes mandatory doubling to min(2*c, 8).
See README.md for the model difference and scripts/generate_chromosome_bound.py
for its regenerated admissible bounds. The inherited PDF describes the parent
model; the ownership and algorithm structure still apply.

Pair-search optimization
------------------------

The capped model now seeds each chromosome's exact pair search with a feasible
history, then searches only for a cheaper one. This changes the search budget,
not the optimal distance or WGD model. Start with ``chromosome_witness`` and
``greedy_witness`` in ``direct_pair/src/pairwise.cpp`` when reviewing this change.
The `optimization note <pair_search_optimization.rst>`_ gives Python-style
pseudocode, the correctness argument and measured results. The inherited PDF
predates this optimization; read this note alongside its pair-search slides.

Ancestral-search optimization
-----------------------------

Sibling joins reuse the pair search's feasible-path helper in
``common/src/feasible_path.hpp``. This supplies a tighter upper bound before
``ancestors/src/nfa_engine.hpp`` runs the exact all-tie search.
Read the `ancestral optimization note <ancestral_search_optimization.rst>`_
for Python-style pseudocode, the offsets/weights accounting and benchmarks.
A subsequent forward-edge filter avoids relaxations that Boost would reject.
Read the `follow-up search note <search_relaxation_optimization.rst>`_ before
reviewing that callback; the backward pass still regenerates all tied labels.
The inherited PDF predates these capped-model search optimizations.

Reducing the search explosion
-----------------------------

Read the `new search review <search_explosion_optimization.rst>`_ before the
native join. OpenFST reduces A followed by capped doubling to a smaller exact
factor. ``ancestors/src/pattern_bound.hpp`` solves a small relaxed sibling
problem to guide Boost A*. ``nfa_engine.hpp`` also indexes suffix transitions by
observed child label. The downward resolver keeps the original primitive order
to preserve historical tie choices. The inherited PDF predates these changes.

Start in Python
---------------

This branch keeps MEDICC's Python package and workflow. The new backend is
medicc.full_fst; compiled sources live under medicc/cext/full_fst, alongside the
native-extension directory already used in medicc2-pre-pub.
This folder is the full-fst-lloh-wgd2-capped-x2-opt branch worktree of medicc2. The original
medicc2 working checkout and the medicc2-pre-pub checkout remain unchanged.

The review route is README.md, this guide, the
`algorithm slides (PDF, 22 slides) <full_fst_algorithms.pdf>`_, then the code.
The slides illustrate the full model, lazy search, retention of all tied
ancestors, and the OpenFST/Boost.Graph boundary using Python-style pseudocode
and a worked example.

Suggested code reading order, all paths relative to this tool folder:

1. medicc/full_fst/core.py: reconstruct() is the public workflow. It encodes
   profiles, obtains a topology from FastME or the caller, reconstructs ancestors,
   and returns Result. This plays the orchestration role of medicc/core.py.
2. medicc/full_fst/scheduler.py: reconstruct() performs ordinary bottom-up joins
   followed by top-down choices, including parallel scheduling.
3. medicc/full_fst/ancestors.py: reconstruct() provides the incremental version
   used by MEDICC NNI. Cache entries own messages shared between neighboring trees.
   This plays the backend role of medicc/ancestors.py. The FULL_TWO sentinel selects
   this backend in the copied MEDICC entry points.
4. medicc/full_fst/model.py and message.py: the small operational vocabulary is
   string(profile), join(left, right), and resolve(parent, child). Read these
   before the native implementation; they expose the mathematical operations.
5. medicc/full_fst/pair_scorer.py: score(x, y) for the shared-ancestor distance,
   or score(x, y, directed=True) for an edge. Both backends share this Python API.
6. medicc/full_fst/_bindings.py: only when you want to inspect the ctypes/C ABI.
   All library paths and ctypes signatures are here.

pipeline.py and integration.py are compatibility imports for existing experiment
scripts. New code uses core.py and ancestors.py. There is only one FULL_TWO token
and one thread-local model/cache namespace, even through the old imports.

The full-FST numerical kernels are C++17. Python owns tree traversal, scheduling,
input/output and NNI orchestration. The new boundary uses ctypes, not Cython;
the inherited FSTlib and acyclic-pruning extensions still use Cython.

Follow one reconstruction
-------------------------

The following is pseudocode; messages remain compressed, not Python word lists::

    for node in bottom_up_order:
        message[node] = model.join(message[left], message[right])
    for node in top_down_order:
        profile[node], cost = model.resolve(profile[parent], message[node])

A leaf starts as model.string(observed_profile). A join searches compatible
primitive histories lazily, finds the best attainable cost, and retains ALL
ancestor words attaining that local cost. OpenFST determinizes that retained
language. The top-down pass selects one optimal word conditional on its parent.

The local all-minimum rule is prune_weight=0. It is not an unpruned global
optimization over every possible internal profile on the tree. Two WGDs are
allowed per directed edge. Equal-cost paths must not be pruned by a >= test.
These are useful invariants to check when reading the C++ search.

Native reading order
--------------------

All paths below start at medicc/cext/full_fst/:

* common/src/automata_types.hpp: the shared integer acceptor representation.
* ancestors/src/nfa_engine.hpp: bottom-up sibling join. Side describes one
  child/history branch; StateKey identifies a product state; SearchState stores
  its current search information; search() explores and retains tight paths.
* ancestors/src/ancestors.cpp: top-down conditional traceback and its C entry
  point. This is the ancestor library's main compilation unit.
* common/src/openfst_language.cpp: conversion to/from OpenFST and its library
  operations. common/src/boost_search.hpp adapts the lazy graph to Boost A*.
* direct_pair/src/pairwise.cpp: the production fixed-profile scorer.

Private implementation headers are included into each engine's compilation
unit. They are not separately linked libraries. The include chain still reflects
the implementation's evolution: shared_engine -> shared_trace -> lazy_engine ->
lazy_trace -> nfa_engine -> ancestors (or direct_trace in the pair engines).
C ABI entry points occur in these included headers as well as the .cpp entry
files. This is intentionally still a reference implementation; splitting every
internal helper into a public header and separate compilation unit is deferred.

Each src/generated/ directory contains large precomputed bound tables. Skip
these for a first algorithm review. They are model constants, not experiment
inputs. budget_2m_core.hpp is a historical filename for primitive operations;
it does not impose a two-million-state search cap. Production search is uncapped.

reference_pair is retained as an independent comparison backend. Its duplicated
native support code is deliberate for validation; it is not a second Python API.
The small pair_engine.hpp support files are retained from the source lineage;
the compiled pair entry point is pairwise.cpp.

Ownership: what the pointers mean
--------------------------------

Model owns the native model handle returned by am_load; close() calls am_free.
Message owns the acceptor returned by am_string, am_import or am_nfa_join;
close() calls am_drop. Both closes are idempotent. A Message retains a Python
reference to its Model, keeping the CDLL available when the acceptor is freed.
PairScorer owns a separate model/cache handle and supports Python's with syntax.

Input messages passed to join/resolve are borrowed, read-only inputs. C++ does
not free them. Python owns the counter arrays and output integers passed by
reference; C++ only writes them during a call. Result strings come from native
thread-local buffers; ctypes c_char_p copies them into Python bytes immediately.
No returned string requires a manual free from Python.

CachedMessage keeps a shared Python reference to a Message. Do not call close()
on a message retained by the NNI cache: other candidate trees may still use it.
A model's mutable search scratch belongs to one worker; immutable messages may
be read by other workers. ExecutionBudget's context manager borrows a worker
permit, and exiting releases that permit. It does not destroy the budget.

A short C++ translation key
--------------------------

* const DFA &child is a borrowed, read-only object, roughly like passing a Python
  object without copying it. The reference must outlive the call using it.
* const DFA *child is a borrowed address in Side. child->len reads its field;
  *child obtains the object. This field does not own or delete the DFA.
* vector<SearchState> nodes is similar to a list of state records. States refer
  to other states by integer IDs. Growing a vector can move its elements, so the
  search reacquires nodes[id] rather than keeping pointers into its old buffer.
* auto asks the compiler to infer a static type; it is not dynamic typing.
* auto &state = nodes[id] refers to the stored record; auto state = nodes[id]
  copies it. Updating the copy does not update the vector.
* A lambda such as [&](int state) { ... } is a local callback. [&] borrows the
  surrounding variables. The callbacks here are used within their owner's scope.
* std::move(value) allows an object to transfer its buffers instead of copying.
  Do not assume the old value still holds its contents afterwards.
* shared_ptr<const vector<int>> shares ownership of an immutable bound table.
  Its last owner frees it. This differs from a borrowed raw pointer.
* RAII means cleanup on scope exit, including exceptions. Stats updates counters
  in its destructor; lock guards release mutexes automatically.
* template<class Nodes> makes a helper work with several statically known node
  containers. Most generic graph machinery is isolated in boost_search.hpp.
* extern "C" exposes stable function names for ctypes. void* is an opaque handle
  at this boundary; Python never dereferences it. Null/error returns become Python
  exceptions in the wrappers.

Which parts are libraries?
--------------------------

OpenFST handles ancestral determinization, weighted union and epsilon removal.
Boost.Graph handles ancestral A* and conditional traceback search. MEDICC keeps
model construction, bounds, bounded lazy composition, tight-path extraction,
acyclic integer minimization/normalization, specialized indexing and the pair
search. Prior measured library alternatives for these retained operations
exceeded the agreed performance budget; the evidence is outside the tool in
experiment 20's library_reuse scripts/results.

What changed for this review
---------------------------

The new core/ancestors names match MEDICC's existing responsibilities. Native
sources moved out of the Python backend into cext. Dynamic loading of Python
files from native folders and two duplicate pair binding implementations were
replaced with an ordinary shared wrapper. Included .cpp implementations became
private .hpp files, large bound tables moved into generated/, and ancestral
search state types gained descriptive names. Algorithms, bounds, tie ordering,
C ABI and compiler options are unchanged.

For validation commands, see README.md. Reproduction drivers, raw timings and
comparison reports belong in experiment 20/review_refactor and the external
medicc_readability_review_20261014 results folder. No datasets are bundled here.

The capped ancestor outlier fallback is described in
``ancestral_outlier_optimization.rst``. It preserves the full prune-zero language
and reuses exact chromosome splitting only after a large search begins.
