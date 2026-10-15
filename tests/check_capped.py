"""High-CN regressions, independent 3-bin oracle and all-tie message checks."""
from pathlib import Path
import itertools
import json
import math
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oracle import get_parts, Primitive, relation, minplus, optimum
from medicc.full_fst.model import Model
from medicc.full_fst.pair_scorer import Model as PairScorer


def main():
    parts, symbols = get_parts()
    scorer = PairScorer()
    reference = PairScorer(backend='reference_pair')
    model = Model()
    checks = dict(targeted=0, pair=0, joins=0, resolves=0)
    try:
        for cn in range(5, 9):
            for x, y in [(f'{cn}X1', '8X2'), (f'{cn}1', '82'), (f'1X{cn}', '2X8')]:
                for pair in [scorer, reference]:
                    assert pair.score(x, y, directed=True) == 1, (x, y)
                message = model.string(y)
                assert model.resolve(x, message) == (y, 1)
                message.close()
                checks['targeted'] += 1
        assert math.isinf(scorer.score('0X5', '1X8', directed=True))
        assert scorer.score('8X1', '8X4', directed=True) == 2
        for template in ['###', '#X##']:
            words = [''.join(next(it) if c == '#' else c for c in template)
                     for digits in itertools.product('012345678', repeat=3) for it in [iter(digits)]]
            matrices = {key: relation(Primitive(parts[key]), words, symbols) for key in ['L_LOH','G','W1']}
            a = minplus(matrices['L_LOH'], matrices['G'])
            t = a
            for _ in range(2):
                t = minplus(minplus(t, matrices['W1']), a)
            k = minplus(t.T, t)
            rng = np.random.default_rng(20261015)
            for _ in range(120):
                i, j = rng.integers(len(words), size=2)
                x, y = words[i], words[j]
                for pair in [scorer, reference]:
                    assert pair.score(x, y) == k[i,j], (x,y,'kernel')
                    assert pair.score(x, y, directed=True) == t[i,j], (x,y,'directed')
                    checks['pair'] += 2
            for _ in range(12):
                i, j, p = rng.integers(len(words), size=3)
                x, y = words[i], words[j]
                left, right = model.string(x), model.string(y)
                joined = model.join(left, right)
                expected = optimum(words, t, {x:0}, {y:0})
                assert joined.language() == expected, (x,y,'all ties')
                checks['joins'] += 1
                # Some random parents have zero at a positive-required bin.
                best = min(t[p,words.index(w)] + cost for w,cost in expected.items())
                if math.isfinite(best):
                    chosen, cost = model.resolve(words[p], joined)
                    assert chosen in expected and cost == best
                    assert t[p,words.index(chosen)] + expected[chosen] == best
                    checks['resolves'] += 1
                joined.close(); left.close(); right.close()
            print(template, checks, flush=True)
    finally:
        scorer.close(); reference.close(); model.close()
    print(json.dumps(dict(status='ok', checks=checks)), flush=True)


if __name__ == '__main__':
    main()
