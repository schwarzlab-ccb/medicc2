from __future__ import annotations

import fstlib
from .acyclic_prune import (
    intersect_prune_acyclic as _intersect_prune_acyclic_native,
    prune_acyclic as _prune_acyclic_native,
)

def intersect_prune_acyclic(
        left: fstlib.Fst,
        right: fstlib.Fst,
        weight: float = 0.0,
        delta: float | None = None,
) -> fstlib.Fst:
    """Intersect two sorted MEDICC acceptors and materialize only pruned paths.

    This experimental fused operation is intentionally narrower than
    ``fstlib.intersect``. Both operands must be standard-arc, epsilon-free,
    input-label-sorted, acyclic, topologically sorted acceptors.
    """
    if not isinstance(left, fstlib.Fst) or not isinstance(right, fstlib.Fst):
        raise TypeError("left and right must be fstlib FSTs")
    if left.arc_type() != "standard" or right.arc_type() != "standard":
        raise ValueError("fused acyclic intersection only supports tropical semiring")
    if delta is None:
        delta = fstlib.DELTA

    native_result = _intersect_prune_acyclic_native(
        left.fst, right.fst, float(weight), float(delta))
    return fstlib.Fst(native_result)

def prune_acyclic(
        ifst: fstlib.Fst,
        weight: float = 0.0,
        delta: float | None = None,
) -> fstlib.Fst:
    """
    Prune a *topologically sorted acyclic standard-arc* FST

     This has the same threshold interpretation as ``fstlib.prune``.  For the
    tropical semiring used by MEDICC2, ``weight=0`` keeps exactly the globally
    minimum-weight successful paths.

    Unlike generic OpenFST pruning, this specialized implementation computes
    distance-to-final values directly by scanning states in descending ID
    order.  It therefore avoids materializing a reversed copy of ``ifst``.

    Parameters
    ----------
    ifst
        An fstlib FST with ``arc_type() == "standard"``.  It must be acyclic,
        topologically sorted, and use the usual dense OpenFST state IDs.
    weight
        Extra pruning threshold in tropical-weight units.
    delta
        OpenFST numerical tolerance.  It is accepted to match fstlib's pruning
        interface; the precomputed topological distances themselves do not
        need iterative convergence.  ``None`` selects fstlib/OpenFST's normal
        default.

    Returns
    -------
    fstlib.MutableFst
        A newly allocated pruned FST.  ``ifst`` is borrowed and never mutated.

    Raises
    ------
    TypeError
        If ``ifst`` is not an fstlib FST.
    ValueError
        If the arc type or graph-shape preconditions are not satisfied.
    RuntimeError
        If the native operation cannot construct an output FST.
    """

    if not isinstance(ifst, fstlib.Fst):
        raise TypeError("ifst must be an fstlib FST.")
    if ifst.arc_type() != "standard":
        raise ValueError("prune_acyclic currently only supports tropical semiring.")

    if delta is None:
        delta = fstlib.DELTA

    native_result = _prune_acyclic_native(ifst.fst, float(weight), float(delta))

    return fstlib.Fst(native_result)



