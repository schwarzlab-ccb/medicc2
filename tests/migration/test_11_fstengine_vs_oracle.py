from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import numpy as np
from medicc.fstengine.model import Model, ModelConfig
with np.load(ROOT / "medicc/test_data/two_bin_oracle.npz") as f:
    words, directed, kernel = f["words"].tolist(), f["directed"], f["undirected"]
with Model() as model:
    for i, x in enumerate(words):
        for j, y in enumerate(words):
            assert model.score(x, y, directed=True) == directed[i, j], (x, y)
            assert model.score(x, y) == kernel[i, j], (x, y)
    a, b = model.string("11"), model.string("22")
    joined = model.join(a, b)
    assert model.resolve("11", joined)[1] == model.resolve("22", joined)[1]  # both ties kept
    for m in (joined, a, b):
        m.close()
with Model(ModelConfig(no_wgd=True)) as model:
    assert model.score("1X1", "2X2", directed=True) == 2
print("PASS step 11: engine == oracle on all two-bin pairs")