"""Independent scores and cache/permutation checks across chromosome boundaries."""
from pathlib import Path
import random
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from medicc.full_fst.pair_scorer import Model

rng = random.Random(20261015)
checks = 0
with Model() as model, Model(backend='reference_pair') as reference:
    cases = [('1X11X111', '2X44X888'), ('8X18X818', '8X28X848'),
             ('0X00X000', '8X88X888'), ('1XX11X', '8XX44X')]
    for _ in range(32):
        cases.append(tuple('X'.join(''.join(rng.choice('012345678') for _ in range(n))
                                   for n in (1, 2, 3)) for _ in range(2)))
    for x, y in cases:
        expected = reference.score(x, y)
        a, b = x.split('X'), y.split('X')
        order = list(range(len(a)))
        rng.shuffle(order)
        px, py = 'X'.join(a[i] for i in order), 'X'.join(b[i] for i in order)
        for left, right in [(x,y), (y,x), (px,py), (py,px), (x,y)]:
            assert model.score(left, right) == expected, (left, right, expected)
            checks += 1
        model.clear()
        assert model.score(x, y) == expected
        checks += 1
print(f'PASS: {len(cases)} independent reference cases; {checks} cache/symmetry/chromosome checks')
