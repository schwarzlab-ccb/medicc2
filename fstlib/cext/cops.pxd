#cython: language_level=3

from libcpp.vector cimport vector

cimport cpywrapfst as fst

cdef extern from "cops.h" nogil:
    cdef void align_std_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, fst.MutableFstClass *)
    cdef fst.WeightClass score_std_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
    cdef fst.WeightClass score_log_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
    cdef fst.WeightClass kernel_score_std_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
    cdef fst.WeightClass kernel_score_log_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
    cdef fst.WeightClass multi_score_std_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
    cdef fst.WeightClass multi_kernel_score_std_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)