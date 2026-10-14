"""Exact pair scores; the two independent native backends share this binding."""

import ctypes as C
import os
from ._bindings import build_directory, load_pair_library


class Model:
    """Own a pair scorer and its cache. Use a separate instance in each thread.

    The reference backend is an independent correctness check. Normal workflows
    use direct_pair. Select a backend explicitly or with MEDICC_FULL_FST_PAIR_BACKEND.
    """

    def __init__(self, bounds="local", max_wgds=2, *, backend=None):
        if bounds != "local" or max_wgds != 2:
            raise ValueError("Fixed full L_LOH two-WGD model")
        self.backend = backend or os.environ.get(
            "MEDICC_FULL_FST_PAIR_BACKEND", "direct_pair"
        )
        if self.backend not in {"direct_pair", "reference_pair"}:
            raise ValueError("Expected direct_pair or reference_pair")
        self.h = None
        self.lib = load_pair_library(self.backend)
        table = build_directory(self.backend) / "ancestral_parts.txt"
        self.h = self.lib.am_load(str(table).encode())
        if not self.h:
            raise RuntimeError(self.lib.am_error().decode())
        # Owned by Python; borrowed by C++ for the duration of score().
        self.counters = (C.c_double * 6)()

    def score(self, x, y, directed=False):
        if not self.h:
            raise ValueError("Closed model")
        if "\0" in x or "\0" in y:
            raise ValueError("NUL in input")
        if any(symbol not in "012345678X" for symbol in x + y):
            raise ValueError("Invalid symbol")
        value = self.lib.am_pair_score(
            self.h, x.encode("ascii"), y.encode("ascii"), directed, self.counters
        )
        if value < 0:
            raise RuntimeError(self.lib.am_error().decode())
        return value

    def clear(self):
        """Drop reusable scores at a row boundary and reset the work counters."""
        if self.h and self.backend == "direct_pair":
            self.lib.am_pair_clear(self.h)
        for i in range(len(self.counters)):
            self.counters[i] = 0

    def stats(self):
        return list(self.counters)

    def close(self):
        """Release the native allocation once; repeated calls are harmless."""
        if self.h:
            self.lib.am_free(self.h)
            self.h = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        if getattr(self, "h", None):
            self.close()
