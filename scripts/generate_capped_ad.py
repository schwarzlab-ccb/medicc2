#!/usr/bin/env python3
"""Generate the exact compact A D factor with the installed OpenFST API.

A is the existing L_LOH G factor; D doubles and saturates at copy number 8.
Output is one AD block that can be appended to ancestral_parts.txt. This is an
asset-generation step: production inference does not run determinization.
"""
import argparse
from collections import deque
import json
from pathlib import Path
import sys


def read_a(path, fst):
    lines = iter(path.read_text().splitlines())
    for header in lines:
        if not header.strip():
            continue
        name, count, start, arc_count = header.split()
        finals = list(map(int, next(lines).split()))
        assert len(finals) == int(count)
        rows = [tuple(map(int, next(lines).split())) for _ in range(int(arc_count))]
        if name != 'A':
            continue
        result = fst.VectorFst()
        for _ in finals:
            result.add_state()
        result.set_start(int(start))
        for state, weight in enumerate(finals):
            result.set_final(state, fst.Weight('tropical', weight))
        for state, source, target, weight, dest in rows:
            result.add_arc(state, fst.Arc(source + 1, target + 1,
                                          fst.Weight('tropical', weight), dest))
        return result
    raise ValueError('Missing A factor')


def count(fst):
    return dict(states=fst.num_states(), arcs=sum(fst.num_arcs(s) for s in fst.states()))


def export(fst):
    # Canonical BFS numbering: deterministic paired labels give a stable traversal
    # order independent of OpenFST's internal minimization state numbering.
    start = fst.start()
    ids = {start: 0}
    queue = deque([start])
    finals, rows = [], []
    while queue:
        state = queue.popleft()
        final = float(fst.final(state))
        assert final >= 0 and final == int(final)
        finals.append(int(final))
        arcs = sorted(fst.arcs(state), key=lambda arc: (arc.ilabel, arc.olabel))
        labels = [(arc.ilabel, arc.olabel) for arc in arcs]
        assert len(labels) == len(set(labels)), 'Expected deterministic encoded labels'
        for arc in arcs:
            assert 1 <= arc.ilabel <= 10 and 1 <= arc.olabel <= 10
            weight = float(arc.weight)
            assert weight >= 0 and weight == int(weight)
            if arc.nextstate not in ids:
                ids[arc.nextstate] = len(ids)
                queue.append(arc.nextstate)
            rows.append((ids[state], arc.ilabel - 1, arc.olabel - 1,
                         int(weight), ids[arc.nextstate]))
    assert len(ids) == fst.num_states()
    # Both search factorization and its lower bounds rely on this reset.
    for state in range(len(finals)):
        assert finals[state] == 0
        assert [r for r in rows if r[0] == state and r[1] == 9] == [(state, 9, 9, 0, 0)]
    return (f'AD {len(finals)} 0 {len(rows)}\n' +
            ' '.join(map(str, finals)) + '\n' +
            ''.join(' '.join(map(str, row)) + '\n' for row in rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True, help='Existing ancestral_parts.txt with factor A')
    ap.add_argument('--output', type=Path, required=True, help='Write one generated AD text block')
    ap.add_argument('--tool', type=Path, help='Tool checkout providing fstlib, if not already importable')
    args = ap.parse_args()
    if args.tool:
        sys.path.insert(0, str(args.tool))
    import fstlib.cext.pywrapfst as fst
    a = read_a(args.source, fst)
    d = fst.VectorFst()
    d.add_state(); d.set_start(0); d.set_final(0)
    for symbol in range(10):
        output = 9 if symbol == 9 else min(2 * symbol, 8)
        d.add_arc(0, fst.Arc(symbol + 1, output + 1, fst.Weight('tropical', 0), 0))
    a.arcsort('olabel'); d.arcsort('ilabel')
    composed = fst.compose(a, d)
    composed.connect()
    mapper = fst.EncodeMapper(arc_type='standard', encode_labels=True, encode_weights=False)
    composed.encode(mapper)
    # Pair encoding is bijective on labels; tropical determinization preserves
    # the minimum weight for each encoded word, including final weights.
    deterministic = fst.determinize(composed)
    compact = deterministic.copy()
    compact.minimize()
    assert fst.equivalent(deterministic, compact), 'Weighted equivalence check failed'
    summary = dict(composed=count(composed), deterministic=count(deterministic),
                   compact=count(compact), weighted_equivalence=True)
    compact.decode(mapper)
    text = export(compact)
    # Check the actual canonical text export too, not only OpenFST's object.
    lines = text.splitlines()
    _, n, start, _ = lines[0].split()
    exported = fst.VectorFst()
    for _ in range(int(n)):
        exported.add_state()
    exported.set_start(int(start))
    for state, weight in enumerate(map(int, lines[1].split())):
        exported.set_final(state, fst.Weight('tropical', weight))
    for row in lines[2:]:
        state, source, target, weight, dest = map(int, row.split())
        exported.add_arc(state, fst.Arc(source + 1, target + 1,
                                        fst.Weight('tropical', weight), dest))
    exported.encode(mapper)
    assert fst.equivalent(deterministic, exported), 'Text export changed the relation'
    summary['export_weighted_equivalence'] = True
    args.output.write_text(text)
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
