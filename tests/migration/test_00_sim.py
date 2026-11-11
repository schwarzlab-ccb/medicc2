from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import medicc.sim
tree = medicc.sim.rcoal(4)
assert sorted(c.name for c in tree.get_terminals()) == ["diploid", "sp_0001", "sp_0002", "sp_0003", "sp_0004"]
print("Pass step 0.2")