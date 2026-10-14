"""Full-FST implementation of MEDICC's incremental NNI interfaces.

Only the explicit FULL_TWO token enables the replacement. Ordinary OpenFST
callers keep the original functions. The copied MEDICC core dispatches directly to these functions.
"""

import math
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

from .model import Model
import medicc
import medicc.ancestors
import medicc.core

FULL_TWO = object()
local = threading.local()


class Genome(str):
    """Immutable profile with the two preparation methods used by nni_mode.

    This is not an OpenFST object. The adapter consumes it as a string. A linear
    string acceptor is already ordered, so copying/sorting requires no change.
    """

    def copy(self):
        return self

    def arcsort(self, order):
        if order != "ilabel":
            raise ValueError(order)
        return self


class CachedMessage:
    """Share one immutable native message between neighbor evaluations.

    `native` is a Python owner of a C++ pointer, not a pointer copy. Retaining
    this wrapper retains that owner and its originating model. Never close a
    shared message manually: its destructor releases C++ memory once the last
    Python reference is gone. Native algorithms only read input messages.
    """

    def __init__(self, native):
        self.native = native

    def arcsort(self, order):
        if order != "ilabel":
            raise ValueError(order)
        return self  # Native DFA transitions are permanently label-indexed.


class Cache(dict):
    """Upper messages plus prior down-pass choices for exact reuse."""

    def __init__(self, values=()):
        super().__init__(values)
        self.profiles = {}
        self.parent_profiles = {}


def model():
    # Every worker owns separate mutable native search scratch and cost cache.
    if not hasattr(local, "model"):
        local.model = Model(state_limit=0, max_wgds=2)
        local.directed = lru_cache(maxsize=None)(
            lambda a, b: local.model.scorer.score(*Model.shorten(a, b), directed=True)
        )
    return local.model


def validate_profiles(profiles, normal_name):
    template = lambda s: "".join("X" if c == "X" else "#" for c in s)
    expected = template(profiles[normal_name])
    if any(
        not isinstance(s, str)
        or template(s) != expected
        or any(c not in "012345678X" for c in s)
        for s in profiles.values()
    ):
        raise ValueError("Invalid reconstructed profile or chromosome grid")


def reconstruct(
    tree,
    samples_dict,
    upper_pass_fst,
    lower_pass_fst,
    normal_name,
    prune_weight=0,
    n_cores=None,
    upper_cache=False,
    old_uppass_cache=None,
    nni_move=None,
    **unused,
):
    if (
        upper_pass_fst is not FULL_TWO
        or lower_pass_fst is not FULL_TWO
        or prune_weight != 0
    ):
        raise ValueError(
            "This adapter requires full two-WGD on both passes and prune(0)"
        )
    if tree.root.name != normal_name or len(tree.root.clades) != 1:
        raise ValueError(
            "Expected MEDICC search shape: one child below the normal root"
        )
    nodes = list(tree.find_clades(order="preorder"))
    names = [n.name for n in nodes]
    if None in names or len(names) != len(set(names)):
        raise ValueError("Node names must be unique")
    internals = [n for n in nodes if n.clades and n is not tree.root]
    if any(len(n.clades) != 2 for n in internals):
        raise ValueError("Expected binary tumour subtree")
    if set(samples_dict) != {normal_name} | {n.name for n in tree.get_terminals()}:
        raise ValueError("Observed sample names do not match seed tips")
    parent = {c.name: n.name for n in nodes for c in n.clades}
    by_name = {n.name: n for n in nodes}
    cache = Cache(old_uppass_cache or {})
    if old_uppass_cache is None:
        dirty = {n.name for n in internals}
    else:
        if nni_move is None:
            raise ValueError("Incremental reconstruction requires an NNI move")
        dirty = set()
        v = nni_move.v_name
        while v != normal_name:
            dirty.add(v)
            v = parent[v]

    def join(name):
        m = model()
        a, b = by_name[name].clades
        # Temporary singleton messages are Python-owned for this call. Cached
        # messages are borrowed read-only; m.join allocates a new output DFA.
        x = (
            m.string(samples_dict[a.name])
            if a.name in samples_dict
            else cache[a.name].native
        )
        y = (
            m.string(samples_dict[b.name])
            if b.name in samples_dict
            else cache[b.name].native
        )
        try:
            return name, CachedMessage(m.join(x, y))
        finally:
            # Only scalar bounds are discarded. No ancestral states or ties.
            m.scorer.clear()

    if old_uppass_cache is None and (n_cores or 1) > 1:
        heights = {}
        levels = {}
        for n in tree.find_clades(order="postorder"):
            heights[n.name] = (
                0 if not n.clades else 1 + max(heights[c.name] for c in n.clades)
            )
            if n.name in dirty:
                levels.setdefault(heights[n.name], []).append(n.name)
        with ThreadPoolExecutor(max_workers=n_cores) as pool:
            for level in sorted(levels):
                cache.update(pool.map(join, levels[level]))
    else:
        for n in reversed(internals):
            if n.name in dirty:
                key, value = join(n.name)
                cache[key] = value

    profiles = {k: Genome(v) for k, v in samples_dict.items()}
    m = model()
    for n in internals:
        p = profiles[parent[n.name]]
        old = old_uppass_cache
        # Reuse a down-pass choice only when BOTH its full upper message and
        # its resolved parent are identical. Otherwise recompute from scratch.
        if (
            isinstance(old, Cache)
            and cache[n.name] is old.get(n.name)
            and old.parent_profiles.get(n.name) == p
        ):
            profiles[n.name] = old.profiles[n.name]
        else:
            profiles[n.name] = Genome(m.resolve(p, cache[n.name].native)[0])
        cache.parent_profiles[n.name] = p
    cache.profiles = profiles
    if unused.get("validate_result", True):
        validate_profiles(profiles, normal_name)
    m.scorer.clear()
    return (
        (profiles, cache) if upper_cache or old_uppass_cache is not None else profiles
    )


def update_lengths(
    tree, fst, ancestor_fsa, normal_name="internal_0", log_progress=True
):
    if fst is not FULL_TWO:
        raise ValueError("Wrong branch-score model")
    model()
    for parent in tree.find_clades():
        for child in parent.clades:
            value = local.directed(ancestor_fsa[parent.name], ancestor_fsa[child.name])
            if not math.isfinite(value):
                raise ValueError("Infeasible reconstructed edge")
            child.branch_length = value
    tree.root.branch_length = 0
    # Scope scalar memoization to one topology; no capacity or search cutoff.
    local.directed.cache_clear()
    local.model.scorer.clear()
