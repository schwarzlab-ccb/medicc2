"""Compatibility imports; new code should use medicc.full_fst.core."""

from .core import (
    Result,
    encode,
    distance_matrix,
    normalize_tree,
    infer_tree,
    decode,
    reconstruct,
)
