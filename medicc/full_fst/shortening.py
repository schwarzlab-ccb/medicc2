"""Precompute sample changes; pair shortening is the union of their boundaries."""

import numpy as np


def prepare(genome):
    values = np.frombuffer(genome.encode("ascii"), dtype=np.uint8)
    keep = np.empty(len(values), dtype=np.bool_)
    if len(values):
        keep[0] = True
        keep[1:] = values[1:] != values[:-1]
        keep[values == ord("X")] = True
    return values, keep


def pair(left, right):
    a, ac = left
    b, bc = right
    if len(a) != len(b):
        raise ValueError("Length mismatch")
    keep = ac | bc
    return a[keep].tobytes().decode("ascii"), b[keep].tobytes().decode("ascii")
