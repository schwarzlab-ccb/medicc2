# cython: language_level=3
# cython: embedsignature=True
# distutils: language=c++

from libcpp.memory cimport shared_ptr
from cython.operator cimport dereference as deref
from fstlib.cext.pywrapfst cimport Fst, MutableFst, _init_MutableFst
cimport cpywrapfst as fst

cdef extern from "acyclic_prune.h" namespace "medicc_fstlib_extension" nogil:
    fst.MutableFstClass* prune_acyclic_std(
        const fst.FstClass& input,
        float weight_threshold,
        float delta,
    ) except +
    fst.MutableFstClass* intersect_prune_acyclic_std(
        const fst.FstClass& left,
        const fst.FstClass& right,
        float weight_threshold,
        float delta,
    ) except +

cpdef MutableFst prune_acyclic(
        Fst ifst,
        float weight=0.0,
        float delta=0.0009765625, # 1 / 1024
):
    cdef shared_ptr[fst.FstClass] input_owner = ifst._fst
    cdef const fst.FstClass* input_ptr = input_owner.get()
    cdef fst.MutableFstClass* output = NULL

    if input_ptr == NULL:
        raise RuntimeError("fstlib Fst contained no native FstClass pointer.")

    with nogil:
        output = prune_acyclic_std(deref(input_ptr), weight, delta)

    if output == NULL:
        raise RuntimeError("native acyclic pruning returned no output FST.")

    return _init_MutableFst(output)

cpdef MutableFst intersect_prune_acyclic(
        Fst left,
        Fst right,
        float weight=0.0,
        float delta=0.0009765625,
):
    """Native fused intersection and exact acyclic pruning."""
    cdef shared_ptr[fst.FstClass] left_owner = left._fst
    cdef shared_ptr[fst.FstClass] right_owner = right._fst
    cdef const fst.FstClass* left_ptr = left_owner.get()
    cdef const fst.FstClass* right_ptr = right_owner.get()
    cdef fst.MutableFstClass* output = NULL

    if left_ptr == NULL or right_ptr == NULL:
        raise RuntimeError("fstlib operand contained no native FstClass pointer.")

    with nogil:
        output = intersect_prune_acyclic_std(
            deref(left_ptr), deref(right_ptr), weight, delta)

    if output == NULL:
        raise RuntimeError("native fused intersection returned no output FST.")

    return _init_MutableFst(output)


