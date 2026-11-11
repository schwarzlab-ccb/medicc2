from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import math
from oracle import Oracle

o = Oracle()
p = o.parts
assert {"L_LOH", "G", "W1"} <= set(p)
# single primitives (model = a chain of one or two FSTs)
assert o.score("5", "8", [p["W1"]]) == 1                 # one capped WGD: 5 -> min(10, 8)
assert o.score("22", "11", [p["L_LOH"]]) == 1            # one interval loss
assert o.score("2X2", "1X1", [p["L_LOH"]]) == 2          # X breaks intervals
A = [p["L_LOH"], p["G"]]
assert math.isinf(o.score("0", "1", A))                  # 0 -> 1: gains cannot resurrect a lost allele
assert o.score("1", "0", A) == 1 and o.score("1", "2", A) == 1
assert o.pair("0", "1", A, A) == 1                       # shared ancestor, zero WGDs: 1 loses one copy toward 0
# the full edge model A (W1 A)^2: at most two WGDs
assert o.directed("1X1", "8X8") == 4                     # 2 gains + 2 WGDs
assert math.isinf(o.directed("0", "1")) and o.kernel("0", "1") == 1
assert o.kernel("0X0", "7X7") == 8 and o.kernel("1X1", "8X8") == o.kernel("8X8", "1X1")
assert o.optimum({"1": 0}, {"2": 0}) == {"1": 1, "2": 1}       # both tied ancestors survive
assert o.optimum({"11": 0}, {"22": 0}) == {"11": 1, "22": 1}   # NOT {"12","21"}: correlations kept
print("PASS step 2")
