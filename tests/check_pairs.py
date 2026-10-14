#!/usr/bin/env python3
"""Exhaustive short-string checks against the independent literal primitive oracle."""
import itertools
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
from oracle import get_parts, Primitive, relation, minplus
import numpy as np

def main():
    parts, symbols = get_parts()
    from medicc.full_fst.pair_scorer import Model
    models = [Model(backend=name) for name in ["direct_pair", "reference_pair"]]
    count = 0
    try:
        for template in ['##', '#X#']:
            words = [''.join(next(it) if c == '#' else c for c in template)
                     for digits in itertools.product('012345678', repeat=2) for it in [iter(digits)]]
            p = {k: relation(Primitive(parts[k]), words, symbols) for k in ['L_LOH','G','W1']}
            a = minplus(p['L_LOH'], p['G']); t = a
            for _ in range(2): t = minplus(minplus(t, p['W1']), a)
            k = minplus(t.T, t)
            for m in models:
                for i, x in enumerate(words):
                    for j, y in enumerate(words):
                        assert m.score(x, y) == k[i, j], (template, x, y, 'kernel')
                        assert m.score(x, y, directed=True) == t[i, j], (template, x, y, 'directed')
                        count += 2
                    m.clear()
        for m in models:
            assert m.score('00014', '13207') == 5
            assert m.score('0X0', '7X7') == 8
            assert m.score('1X1', '8X8', directed=True) == 4
        print(f'PASS: {count} oracle comparisons across both pair backends')
    finally:
        for m in models: m.close()
if __name__ == '__main__': main()
