from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import pandas as pd
from medicc.fstengine import ModelConfig, pairwise_distances, reconstruct
samples = {"b": "22", "a": "11"}
one = pairwise_distances(samples, ModelConfig(), n_cores=1)
two = pairwise_distances(samples, ModelConfig(), n_cores=2)
pd.testing.assert_frame_equal(one, two)
assert list(one.index) == ["b", "a"] and one.loc["b", "a"] == 1
results = [reconstruct({"r": ["a", "b"]}, samples, "r", "11", ModelConfig(), n_cores=n) for n in (1, 2)]
assert results[0] == results[1] == {"b": "22", "a": "11", "r": "11"}
print("PASS step 14")