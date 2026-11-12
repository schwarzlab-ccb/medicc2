"""Freeze the oracle's answers on every two-bin profile.

Writes medicc/test_data/two_bin_oracle.npz with
    words       the 81 words "00" .. "88"
    directed    81x81, directed[i, j] = T(words[i], words[j])   (parent -> child)
    undirected  81x81, undirected[i, j] = K(words[i], words[j]) (pair distance)
inf means "impossible". Run once:  python -B tests/make_golden.py
"""
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
from oracle import MAX_CN, Oracle, _total, _words

OUTPUT = ROOT / "medicc/test_data/two_bin_oracle.npz"


def main():
    start = time.time()
    o = Oracle()
    words = [f"{a}{b}" for a in range(MAX_CN + 1) for b in range(MAX_CN + 1)]
    index = {w: i for i, w in enumerate(words)}
    n = len(words)

    # T: push each parent once through the model and read off every child with its cost.
    directed = np.full((n, n), np.inf)
    for x in words:
        for y, cost in _words(o.children(o.language({x: 0})), o.sym).items():
            directed[index[x], index[y]] = cost

    # K: the ancestor language of each word once, then one intersection per pair.
    ancestors = [o.ancestors(o.language({w: 0})) for w in words]
    undirected = np.full((n, n), np.inf)
    for i in range(n):
        for j in range(i, n):
            undirected[i, j] = undirected[j, i] = _total(ancestors[i] * ancestors[j])

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUTPUT, words=np.array(words), directed=directed, undirected=undirected)
    print(f"wrote {OUTPUT.relative_to(ROOT)} in {time.time() - start:.0f}s")


if __name__ == "__main__":
    main()
