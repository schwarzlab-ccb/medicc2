"""Full L_LOH two-WGD experimental API: string -> message, join two messages, resolve one ancestor.

Read core.reconstruct for the public workflow and ancestors.reconstruct for the
bottom-up / top-down traversal. ABI declarations live in _bindings.py.
"""

import ctypes as C
import json
from contextlib import nullcontext
from ._bindings import build_directory, load_ancestor_library, load_budget_library
from .message import Message, INF
from .pair_scorer import Model as PairScorer

BUILD = build_directory("ancestors")


class Model:
    """Own one native model; use a different instance for each concurrent thread."""

    def __init__(self, state_limit=0, max_wgds=2):
        if state_limit != 0 or max_wgds != 2:
            raise ValueError(
                "This L_LOH build is uncapped and permits at most two WGDs"
            )
        self.execution_budget = None
        self.handle = None
        self.scorer = None
        self.lib = load_ancestor_library()
        self.handle = self.lib.am_load(str(BUILD / "ancestral_parts.txt").encode())
        if not self.handle:
            raise RuntimeError(self.lib.am_error().decode())
        self.scorer = PairScorer(bounds="local", max_wgds=2)
        # Python owns these fixed arrays. C++ borrows their addresses during calls.
        self.counters = (C.c_double * 5)()
        self.down_counters = (C.c_double * 3)()

    def close(self):
        """Idempotent: release each native allocation once, even after an error."""
        if self.scorer is not None:
            self.scorer.close()
        if self.handle:
            self.lib.am_free(self.handle)
            self.handle = None

    def __del__(self):
        if getattr(self, "handle", None):
            self.close()

    @staticmethod
    def shorten(a, b):
        """Collapse equal adjacent column pairs only for the witness pair scorer.

        Reconstruction messages retain every original bin. X boundaries are kept.
        """
        if len(a) != len(b):
            raise ValueError("Length mismatch")
        keep = [
            i
            for i in range(len(a))
            if i == 0
            or a[i] == "X"
            or a[i - 1] == "X"
            or (a[i], b[i]) != (a[i - 1], b[i - 1])
        ]
        return "".join(a[i] for i in keep), "".join(b[i] for i in keep)

    def string(self, genome):
        if not self.handle:
            raise ValueError("Closed model")
        if not isinstance(genome, str) or any(c not in "012345678X" for c in genome):
            raise ValueError("Expected digits 0..8 and X")
        # C++ returns a newly allocated DFA. Message becomes its sole owner.
        return Message(self, self.lib.am_string(genome.encode()))

    def join(self, left, right):
        with self.execution_budget or nullcontext():
            return self._join(left, right)

    def _join(self, left, right):
        """Return ALL minimum-cost ancestors for these two retained messages."""
        if not self.handle or not left.handle or not right.handle:
            raise ValueError("Closed handle")
        if left.offset + right.offset >= INF:
            raise ValueError("Combined cost exceeds supported integer range")
        x, wx = left.choose()
        y, wy = right.choose()
        # The chosen child profiles are feasible. Their exact kernel score plus
        # stored subtree costs is an upper bound, never an empirical cutoff.
        upper = self.scorer.score(*self.shorten(x, y)) + wx + wy
        upper = min(upper, INF)
        handle = self.lib.am_nfa_join(
            self.handle, left.handle, right.handle, int(upper), 0, 0, self.counters
        )
        return Message(self, handle)

    def resolve(self, parent, child):
        with self.execution_budget or nullcontext():
            return self._resolve(parent, child)

    def _resolve(self, parent, child):
        """Choose one optimal child profile conditional on a fixed parent string."""
        if not self.handle or not child.handle:
            raise ValueError("Closed handle")
        if any(c not in "012345678X" for c in parent):
            raise ValueError("Invalid parent symbol")
        chosen, offset = child.choose()
        upper = self.scorer.score(*self.shorten(parent, chosen), directed=True) + offset
        upper = min(upper, INF)
        cost = C.c_int()
        # byref(cost) means 'let C++ fill this integer'; Python keeps it alive.
        word = self.lib.am_direct_resolve(
            self.handle,
            parent.encode(),
            child.handle,
            int(upper),
            C.byref(cost),
            self.down_counters,
        )
        if word is None:
            raise RuntimeError(self.lib.am_error().decode())
        return word.decode(), cost.value

    def stats(self):
        names = [
            "nfa_expansions",
            "max_nfa_product_states",
            "generated_nfa_transitions",
            "max_retained_nfa_states",
            "max_output_subset_states",
        ]
        result = dict(zip(names, self.counters))
        result.update(
            zip(
                ["down_expansions", "down_peak_states", "down_generated"],
                self.down_counters,
            )
        )
        trace = self.lib.am_last_join_trace()
        if trace:
            result["join_trace"] = json.loads(trace)
        return result

    def from_fst(self, fsa, length):
        """Import a deterministic, epsilon-free, fixed-length weighted acceptor.

        Refuse unsupported weights/labels rather than silently changing semantics.
        Caller owns determinization if the existing message is nondeterministic.
        """
        import math

        if fsa.arc_type() != "standard":
            raise ValueError("Expected tropical acceptor")
        f = fsa.fst
        n = f.num_states()
        rows = []

        def weight(w):
            x = float(w)
            if x == math.inf:
                return INF
            if not math.isfinite(x) or x < 0 or x != int(x) or x >= INF:
                raise ValueError("Expected finite nonnegative integer weight")
            return int(x)

        for s in range(n):
            row = [weight(f.final(s))] + [-1, 0] * 10
            seen = set()
            for a in f.arcs(s):
                if not a.ilabel or not a.olabel:
                    raise ValueError("Epsilon arcs unsupported")
                x = fsa.input_symbols().find(a.ilabel)
                y = fsa.output_symbols().find(a.olabel)
                if x != y or x not in list("012345678X"):
                    raise ValueError("Expected copy-number acceptor")
                c = 9 if x == "X" else int(x)
                if c in seen:
                    raise ValueError("Expected deterministic acceptor")
                seen.add(c)
                row[1 + 2 * c] = a.nextstate
                row[2 + 2 * c] = weight(a.weight)
            rows.extend(row)
        finite_finals = [rows[21 * i] for i in range(n) if rows[21 * i] < INF]
        max_arc = max(
            (rows[21 * i + 2 + 2 * c] for i in range(n) for c in range(10)), default=0
        )
        if length * max_arc + max(finite_finals, default=0) >= INF:
            raise ValueError("Acceptor exceeds supported total-cost range")
        data = (C.c_int * len(rows))(*rows)
        return Message(self, self.lib.am_import(length, n, f.start(), 0, data))


class ExecutionBudget:
    """Own the native-work budget shared by tree workers and WGD helpers.

    Entering/exiting borrows/releases a worker permit; it does not create or
    destroy the budget. The scheduler keeps this owner alive until workers join.
    """

    def __init__(self, workers):
        self.lib = load_budget_library()
        self.handle = self.lib.am_budget_new(workers)
        if not self.handle:
            raise RuntimeError(self.lib.am_error().decode())

    def __enter__(self):
        if not self.lib.am_worker_enter(self.handle):
            raise RuntimeError(self.lib.am_error().decode())
        return self

    def __exit__(self, *exc):
        self.lib.am_worker_leave()

    @property
    def peak(self):
        return self.lib.am_budget_peak(self.handle)

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.am_budget_free(self.handle)
