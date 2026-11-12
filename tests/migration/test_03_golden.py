from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np

with np.load(ROOT / "medicc/test_data/two_bin_oracle.npz") as f:
    words = f["words"].tolist()
    directed, kernel = f["directed"], f["undirected"]
assert set(words) == {f"{a}{b}" for a in range(9) for b in range(9)}
assert directed.shape == kernel.shape == (81, 81)
assert np.array_equal(kernel, kernel.T)                          # K is symmetric
assert np.all(np.diag(directed) == 0)                            # no change costs 0
assert np.isinf(directed[words.index("00"), words.index("11")])  # cannot gain from 0
assert np.all(kernel <= directed)                                # the child itself can be the ancestor
# K(x, y) = min_z T(z, x) + T(z, y); for two-bin words every ancestor z is a two-bin word too
assert np.array_equal(kernel, np.min(directed[:, :, None] + directed[:, None, :], axis=0))
print("PASS step 3")
