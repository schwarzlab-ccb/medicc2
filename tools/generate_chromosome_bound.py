#!/usr/bin/env python3
"""Generate medicc/native/src/chromosome_constants.hpp.

The tables are boundary relaxations of the native model A (D A)^k, k = 0, 1, 2,
with A = L_LOH o G and the capped doubling D(c) = min(2c, 8). Each interval
event contributes an opening and a closing boundary. Compatibility between
boundaries is relaxed, each boundary is minimised independently, and the sum
divided by two (rounded up) is an admissible lower bound on the event count.
L_LOH permits interval loss depth beyond the CN when the output reaches zero;
gain is its non-resurrecting counterpart. WGD costs are added outside the tables.

usage: python tools/generate_chromosome_bound.py
"""
import itertools
from pathlib import Path

import numpy as np

OUTPUT = Path(__file__).resolve().parents[1] / "medicc" / "native" / "src" / "chromosome_constants.hpp"
INF = 10000  # "no path": a big int, not float infinity, because the tables are int32


# Matrix product in the tropical semiring: * becomes +, sum becomes min.
#   out[i, j] = min over k of a[i, k] + b[k, j]
# If a is "cost from state i to state k" and b is "cost from k to j", out is "cost from i to j
# through any middle state k": the two steps chained into one (like composing two FSTs, costs only).
def minplus(a, b):
    out = np.full((len(a), b.shape[1]), INF, dtype=np.int32)
    # One middle state k at a time: a[:, k:k+1] + b[k:k+1, :] is the full i x j table of
    # a[i, k] + b[k, j]; keep the smallest value seen so far.
    for k in range(a.shape[1]):
        np.minimum(out, a[:, k:k+1] + b[k:k+1, :], out=out)
    return out


# The interval functions answer, for one column: how many events of this FST must cover it to
# turn copy number `source` into `target`? They return (lo, hi), the fewest and the most events,
# or None if no number of events can do it. The number of events that cover a column is its
# "depth". The ranges must match the model, not the biology: in the model 8 is exactly 8.
# A hi that is too big is safe (the bound only gets weaker); a hi that is too small is not.

def interval_l_loh(source, target):
    """Losses (L_LOH): each event takes away 1. 0 is a floor: extra losses keep it at 0 (LOH)."""
    if target > source:
        return None  # a loss never goes up
    # 3 -> 0: at least 3 losses, at most 8 (L_LOH has 8 steps, max_pre_wgd_losses in factory.py).
    # 3 -> 1: exactly 2 losses.
    return (source, 8) if target == 0 else (source-target, source-target)


def interval_g(source, target):
    """Gains (G): each event adds 1. A 0 is never gained, and a gain never makes a 0."""
    if source == target == 0:
        return 0, 8  # a gain over a 0 does nothing: 0 to 8 gains (G has 8 steps) can cover it
    if 0 < source <= target:
        return target-source, target-source  # 1 -> 3: exactly 2 gains
    return None  # 0 -> 2 (no gain from 0), 2 -> 0 (a gain never makes a 0), 3 -> 1 (going down)


def boundary_cost(caps):
    """The fewest event boundaries for the (lo, hi) of each column (one column: a chromosome end;
    two columns: neighbours). INF if a column is impossible."""
    if any(c is None for c in caps):
        return INF
    if len(caps) == 1:
        # A chromosome end: every event that covers this column must start (or end) here.
        return caps[0][0]
    # Two neighbours: an event that covers one but not the other starts or ends between them.
    # The depths can be chosen inside their ranges, so only the gap between the ranges counts:
    #   (2,2) and (2,2) -> 0: the same 2 events cover both columns.
    #   (2,2) and (0,0) -> 2: 2 events end between them.
    #   (3,8) and (2,2) -> 1: the left needs >= 3, the right allows 2, so >= 1 ends between them.
    (lo, hi), (lo2, hi2) = caps
    return max(0, lo-hi2, lo2-hi)


# The primitive tables: one per FST, out[i, j] = the fewest boundaries to turn state i into
# state j with that FST alone (INF if it cannot). A state is one copy number (a chromosome end)
# or a pair (left, right) of neighbouring copy numbers; each column of the state is checked
# with the FST's interval function, then boundary_cost combines them.

def primitive_l_loh(states):
    out = np.full((len(states), len(states)), INF, dtype=np.int32)
    for i, source in enumerate(states):
        for j, target in enumerate(states):
            out[i, j] = boundary_cost([interval_l_loh(s, t) for s, t in zip(source, target)])
    return out


def primitive_g(states):
    out = np.full((len(states), len(states)), INF, dtype=np.int32)
    for i, source in enumerate(states):
        for j, target in enumerate(states):
            out[i, j] = boundary_cost([interval_g(s, t) for s, t in zip(source, target)])
    return out


# The tables for the whole model. endpoint=True: states are one copy number 0..8 (9 states);
# endpoint=False: states are pairs (left, right) of neighbours (81 states). capped=True: the
# doubling is D(c) = min(2c, 8); capped=False: no doubling above 4 (main only uses capped=True).
def tables(endpoint, capped):
    if endpoint:
        columns = 1  # one column: [(0,), (1,), ..., (8,)], 9 states
    else:
        columns = 2  # two neighbours: [(0,0), (0,1), ..., (8,8)], 81 states
    states = list(itertools.product(range(9), repeat=columns))
    ids = {s: i for i, s in enumerate(states)}
    # A = L_LOH o G: first the losses, then the gains.
    a = minplus(primitive_l_loh(states), primitive_g(states))
    # D: one WGD doubles every column, at no cost here (WGD costs are added outside the tables).
    d = np.full_like(a, INF)
    for i, source in enumerate(states):
        if capped or max(source) <= 4:
            target = tuple(min(8, 2*c) for c in source)
            d[i, ids[target]] = 0
    # One branch (parent -> child) with k WGDs: A, ADA, ADADA for k = 0, 1, 2.
    branches = [a]
    for _ in range(2):
        branches.append(minplus(minplus(branches[-1], d), a))
    # Two branches from a common ancestor z, one with ka WGDs to x and one with kb WGDs to y:
    #   minplus(bx.T, by)[x, y] = min over z of bx[z, x] + by[z, y]
    # 9 tables, index 3*ka + kb.
    return np.stack([minplus(x.T, y) for x in branches for y in branches])


# Write the tables as C++ constants: END[9][9*9] (one copy number) and EDGE[9][81*81]
# (neighbour pairs); each row is one (ka, kb) table, flattened.
def main():
    text = '// Generated by tools/generate_chromosome_bound.py.\n'
    text += '// Capped D(c)=min(2*c,8); branch index 3*ka+kb. WGD costs excluded.\n'
    text += '#pragma once\nnamespace chromosome_bound {\n'
    for name, endpoint in [('END', True), ('EDGE', False)]:
        table = tables(endpoint, True).reshape(9, -1)
        assert 0 <= table.min() <= table.max() < 255
        text += f'static constexpr unsigned char {name}[9][{table.shape[1]}] = {{\n'
        text += ''.join(' {' + ','.join(map(str, row)) + '},\n' for row in table)
        text += '};\n'
    text += '\n} // namespace chromosome_bound\n'
    OUTPUT.write_text(text)
    print(f'Wrote {OUTPUT}')


if __name__ == '__main__':
    main()
