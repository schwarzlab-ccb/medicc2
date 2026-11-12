"""Build the fstengine's model parts with medicc.factory and write them as text.

Parts:
    L_LOH   losses incl. LOH (medicc.factory)
    G       gains (medicc.factory)
    A       L_LOH o G, made small
    AD      A o D, made small; D = capped doubling c -> min(2c, 8), X -> X, cost 0
            (the engine charges the WGD cost once per genome, not per bin)

File format, one block per part:
    NAME n_states start n_arcs
    final weight of state 0 .. n-1 (one line; NOT_FINAL = 100000000 means "not final")
    source ilabel olabel cost target (one line per arc)
Labels: copy number 0..8, X -> 9 (the engine adds 1 so label 0 stays epsilon)
States are numbered breadth-first from the start state (start = 0), arcs sorted by labels,
so the output is the same on every run.

Usage: python tools/export_model_parts.py --output parts.txt
"""
import argparse
from collections import deque
from pathlib import Path

import fstlib
import medicc

MAX_CN = 8
SEPARATOR = 'X'
NOT_FINAL = 100000000

def doubling(sym):
    """D: every symbol c -> min(2c, MAX_CN), X -> X, cost 0. One state, start and final."""
    d = fstlib.Fst()
    d.set_input_symbols(sym)
    d.set_output_symbols(sym)
    q = d.add_state()
    d.set_start(q)
    d.set_final(q, 0)
    for c in [str(i) for i in range(MAX_CN + 1)] + [SEPARATOR]:
        d.add_arc(q, (c, c if c == SEPARATOR else str(min(2 * int(c), MAX_CN)), 0, q))
    return d

def build_parts():
    sym = medicc.create_symbol_table(MAX_CN, SEPARATOR)
    p = medicc.create_copynumber_fst(symbol_table=sym, sep=SEPARATOR, enable_wgd=False, output_all=True)
    a = fstlib.encode_determinize_minimize(p['L_LOH'] * p['G'])
    ad = fstlib.encode_determinize_minimize(a * doubling(sym))
    return {"L_LOH": p['L_LOH'], 'G': p['G'], 'A': a, 'AD': ad}, sym

def export(name, f, sym):
    """One text block. Breadth-first state numbering from the start state
    File format, one block per part:
        NAME n_states start n_arcs
        final weight of state 0 .. n-1 (one line; NOT_FINAL = 100000000 means "not final")
        source ilabel olabel cost target (one line per arc)
    """
    label = {str(c): c for c in range(MAX_CN + 1)}
    label[SEPARATOR] = MAX_CN + 1
    ids = {f.start(): 0}
    queue = deque([f.start()])
    finals, rows = [], []
    while queue:
        state = queue.popleft()
        final = float(f.final(state))
        if final == float("inf"):
            finals.append(NOT_FINAL)
        else:
            assert final >= 0 and final == int(final), (name, state, final)
            finals.append(int(final))
        arcs = sorted(f.arcs(state), key=lambda arc: (arc.ilabel, arc.olabel, float(arc.weight)))
        for arc in arcs:
            weight = float(arc.weight)
            assert weight >= 0 and weight == int(weight), (name, state, weight)
            if arc.nextstate not in ids:
                ids[arc.nextstate] = len(ids)
                queue.append(arc.nextstate)
            rows.append((ids[state], label[sym.find(arc.ilabel)], label[sym.find(arc.olabel)],
                         int(weight), ids[arc.nextstate]))
    assert len(ids) == f.num_states(), f'{name}: unreachable states'
    return (f'{name} {len(finals)} 0 {len(rows)}\n' + ' '.join(map(str, finals)) + '\n'
            + ''.join(' '.join(map(str, r)) + '\n' for r in rows))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', type=Path, required=True, help='text file to write')
    args = ap.parse_args()
    parts, sym = build_parts()
    args.output.write_text(''.join(export(name, f, sym) for name, f in parts.items()))
    for name, f in parts.items():
        print(f'{name}: {f.num_states()} states, {sum(1 for s in f.states() for _ in f.arcs(s))} arcs')
    print(f'wrote {args.output}')

if __name__ == '__main__':
    main()
