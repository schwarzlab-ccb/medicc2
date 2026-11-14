"""REFERENCE for migration Step 28 (not part of REF's code). Goes to tests/migration/test_28_shortening.py.
The plan's test, plus: the scheduler's distances equal the scores of the full profiles.
"""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
from medicc.fstengine import ModelConfig
from medicc.fstengine.scheduler import pairwise_distances
from medicc.fstengine.shortening import prepare, pair
from medicc.fstengine.model import Model
x, y = "111X22", "122X44"
assert pair(prepare(x), prepare(y)) == ("11X2", "12X4")
assert pair(prepare("1XX1"), prepare("2XX2")) == ("1XX1", "2XX2")
profiles = {"a": "1112X2333", "b": "1222X2211", "c": "2222X4444"}
distances = pairwise_distances(profiles, ModelConfig())
with Model() as model:
    assert model.score(x, y) == model.score(*pair(prepare(x), prepare(y)))
    for a in profiles:
        for b in profiles:
            if a != b:
                assert distances.loc[a, b] == model.score(profiles[a], profiles[b]), (a, b)
print("PASS step 28")
