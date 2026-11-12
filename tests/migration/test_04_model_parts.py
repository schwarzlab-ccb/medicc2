from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

tokens = iter((ROOT / "medicc/objects/ancestral_parts.txt").read_text().split())
seen = set()
for name in tokens:
    seen.add(name)
    states, start, count = (int(next(tokens)) for _ in range(3))
    assert states > 0 and 0 <= start < states
    finals = [float(next(tokens)) for _ in range(states)]
    for _ in range(count):
        s, i, o, w, t = (int(next(tokens)) for _ in range(5))
        assert 0 <= s < states and 0 <= t < states and 0 <= i <= 9 and 0 <= o <= 9 and w >= 0
assert {"L_LOH", "G", "A", "AD"} <= seen
print("PASS step 4")