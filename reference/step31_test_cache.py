"""REFERENCE for migration Step 31 (not part of REF's code). Goes to tests/migration/test_31_cache.py.
The plan's test, plus: a long-lived cache gives the same answers as a fresh one.
"""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import itertools
from medicc.fstengine.model import Model
with Model() as model:
    first = model.score("1221X3", "2112X1")
    model._finish_row()
    assert model.score("1221X3", "2112X1") == first
    # 5120 pairs through one Model (its cache only grows, until the limit clears it) must give
    # the same answers as a fresh cache for every pair.
    words = [''.join(w) for w in itertools.product('1234', repeat=4)]   # 256 words
    expected = {}
    with Model() as fresh:
        for x in words[:20]:
            for y in words:
                expected[x, y] = fresh.score(x, y)
                fresh.clear()   # a fresh cache for every pair: the reference answer
    for x, y in expected:
        assert model.score(x, y) == expected[x, y], (x, y)
print("PASS step 31")
