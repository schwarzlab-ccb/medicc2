"""Column shortening for pair distances.

A column is dropped if it repeats the previous column in both profiles (as in
``medicc.core.shorten_cn_strings``) and is not a chromosome separator.
Per-sample change masks are computed once; a pair keeps the union of both masks.
"""
import numpy as np


def prepare(profile):
    """Return the profile as bytes and the mask of columns that differ from their predecessor."""
    values = np.frombuffer(profile.encode('ascii'), dtype=np.uint8)
    keep = np.empty(len(values), dtype=np.bool_)
    if len(values):
        keep[0] = True
        keep[1:] = values[1:] != values[:-1]
        keep[values == ord('X')] = True
    return values, keep


def pair(left, right):
    """Shorten two prepared profiles of equal length to the columns that either one changes."""
    a, a_keep = left
    b, b_keep = right
    if len(a) != len(b):
        raise ValueError("Profiles have different lengths")
    keep = a_keep | b_keep
    return a[keep].tobytes().decode('ascii'), b[keep].tobytes().decode('ascii')
