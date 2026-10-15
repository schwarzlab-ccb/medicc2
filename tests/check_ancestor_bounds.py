"""Compare complete weighted ancestral languages against a literal FST oracle."""
import ctypes as C
import itertools
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'tests')]
from oracle import get_parts, Primitive, relation, minplus, optimum
from medicc.full_fst.model import Model, ExecutionBudget
from medicc.full_fst.message import Message, INF


def weighted(model, language, offset):
    """Import a trie with costs on arcs AND finals, plus a nonzero offset."""
    nodes = [dict(final=INF, arcs={})]
    for word, value in language.items():
        state = 0
        for char in word:
            symbol = 9 if char == 'X' else int(char)
            if symbol not in nodes[state]['arcs']:
                nodes[state]['arcs'][symbol] = len(nodes)
                nodes.append(dict(final=INF, arcs={}))
            state = nodes[state]['arcs'][symbol]
        nodes[state]['final'] = value
    rows = []
    for node in nodes:
        row = [node['final']] + [-1, 0]*10
        for symbol, dest in node['arcs'].items():
            row[1+2*symbol:3+2*symbol] = [dest, 1]
        rows.extend(row)
    length = len(next(iter(language)))
    data = (C.c_int*len(rows))(*rows)
    message = Message(model, model.lib.am_import(length, len(nodes), 0, offset, data))
    expected = {word: cost+length+offset for word, cost in language.items()}
    assert message.language() == expected
    return message, expected


def main():
    parts, symbols = get_parts()
    rng = np.random.default_rng(20261015)
    model = Model()
    checked = 0
    try:
        for template in ['##', '#X#']:
            words = [''.join(next(it) if c == '#' else c for c in template)
                     for digits in itertools.product('012345678', repeat=2)
                     for it in [iter(digits)]]
            matrices = {key: relation(Primitive(parts[key]), words, symbols)
                        for key in ['L_LOH', 'G', 'W1']}
            a = minplus(matrices['L_LOH'], matrices['G'])
            t = a
            for _ in range(2):
                t = minplus(minplus(t, matrices['W1']), a)
            for case in range(48):
                languages = [{words[int(i)]: int(rng.integers(6))
                              for i in rng.choice(len(words), size=5, replace=False)}
                             for _ in range(2)]
                left, lv = weighted(model, languages[0], 3)
                right, rv = weighted(model, languages[1], 7)
                expected = optimum(words, t, lv, rv)
                try:
                    for workers in [1, 2]:
                        model.execution_budget = ExecutionBudget(workers)
                        joined = model.join(left, right)
                        try:
                            assert joined.language() == expected, (template, case, workers)
                            checked += 1
                        finally:
                            joined.close()
                finally:
                    model.execution_budget = None
                    left.close()
                    right.close()
    finally:
        model.close()
    print(json.dumps(dict(status='ok', weighted_all_tie_joins=checked)), flush=True)


if __name__ == '__main__':
    main()
