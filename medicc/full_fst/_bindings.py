"""Private ctypes boundary: library paths, signatures, and model checks.

Python wrappers own opaque handles; only C++ dereferences them. Declaring pointer
return types here is essential: ctypes otherwise assumes a 32-bit integer.
C++ string results use c_char_p, so ctypes copies them into Python bytes before
another native call can overwrite the thread-local result buffer.
"""

import ctypes as C
from pathlib import Path

NATIVE = Path(__file__).resolve().parents[1] / "cext/full_fst"


def build_directory(backend):
    return NATIVE / backend / "build"


def _load(backend, filename, signatures):
    library = C.CDLL(str(build_directory(backend) / filename))
    for name, arguments, result in signatures:
        function = getattr(library, name)
        function.argtypes = arguments
        function.restype = result
    return library


def _check_model(library):
    if (
        library.am_max_wgds() != 2
        or library.am_lloh_model() != 1
        or library.am_fused_model() != 1
    ):
        raise ValueError("Loaded library is not the full L_LOH two-WGD model")


def load_ancestor_library():
    signatures = [
        ("am_error", [], C.c_char_p),
        ("am_last_join_trace", [], C.c_char_p),
        ("am_load", [C.c_char_p], C.c_void_p),
        ("am_free", [C.c_void_p], None),
        ("am_drop", [C.c_void_p], None),
        ("am_string", [C.c_char_p], C.c_void_p),
        (
            "am_import",
            [C.c_int, C.c_int, C.c_int, C.c_int, C.POINTER(C.c_int)],
            C.c_void_p,
        ),
        ("am_size", [C.c_void_p], C.c_int),
        ("am_start", [C.c_void_p], C.c_int),
        ("am_offset", [C.c_void_p], C.c_int),
        ("am_final", [C.c_void_p, C.c_int], C.c_int),
        ("am_arc", [C.c_void_p, C.c_int, C.c_int, C.POINTER(C.c_int)], C.c_int),
        ("am_choose_profile", [C.c_void_p, C.POINTER(C.c_int)], C.c_char_p),
        (
            "am_nfa_join",
            [
                C.c_void_p,
                C.c_void_p,
                C.c_void_p,
                C.c_int,
                C.c_int,
                C.c_int,
                C.POINTER(C.c_double),
            ],
            C.c_void_p,
        ),
        (
            "am_direct_resolve",
            [
                C.c_void_p,
                C.c_char_p,
                C.c_void_p,
                C.c_int,
                C.POINTER(C.c_int),
                C.POINTER(C.c_double),
            ],
            C.c_char_p,
        ),
        ("am_max_wgds", [], C.c_int),
        ("am_lloh_model", [], C.c_int),
        ("am_fused_model", [], C.c_int),
    ]
    library = _load("ancestors", "libdirect_wgd2.so", signatures)
    _check_model(library)
    return library


def load_budget_library():
    signatures = [
        ("am_budget_new", [C.c_int], C.c_void_p),
        ("am_budget_free", [C.c_void_p], None),
        ("am_worker_enter", [C.c_void_p], C.c_int),
        ("am_worker_leave", [], None),
        ("am_budget_peak", [C.c_void_p], C.c_int),
        ("am_error", [], C.c_char_p),
    ]
    return _load("ancestors", "libdirect_wgd2.so", signatures)


def load_pair_library(backend):
    signatures = [
        ("am_load", [C.c_char_p], C.c_void_p),
        ("am_free", [C.c_void_p], None),
        ("am_error", [], C.c_char_p),
        (
            "am_pair_score",
            [C.c_void_p, C.c_char_p, C.c_char_p, C.c_int, C.POINTER(C.c_double)],
            C.c_double,
        ),
        ("am_max_wgds", [], C.c_int),
        ("am_lloh_model", [], C.c_int),
        ("am_fused_model", [], C.c_int),
        ("am_pair_score_only", [], C.c_int),
        ("am_pair_search_uncapped", [], C.c_int),
    ]
    if backend == "direct_pair":
        signatures.extend(
            [
                ("am_pair_clear", [C.c_void_p], None),
                ("am_pair_cache_size", [C.c_void_p], C.c_size_t),
            ]
        )
    library = _load(backend, "libpair_direct.so", signatures)
    _check_model(library)
    if not (library.am_pair_score_only() and library.am_pair_search_uncapped()):
        raise ValueError("Loaded library is not an uncapped pair scorer")
    return library
