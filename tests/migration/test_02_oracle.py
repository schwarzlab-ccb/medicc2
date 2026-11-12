from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import math
import medicc
from oracle import Oracle, _words

o = Oracle()
p = o.parts
assert {"L_LOH", "G", "D"} <= set(p)
assert len(o.model) == 3                                     # A, ADA, ADADA
# single primitives (model = one chain of one or two FSTs)
assert o.score("5", "8", [[p["D"]]]) == 1                    # one capped WGD: 5 -> min(10, 8)
assert math.isinf(o.score("5", "5", [[p["D"]]]))             # D is mandatory: no identity path
assert o.score("0X3", "0X6", [[p["D"]]]) == 1                # cost 1 once per genome, X kept
assert o.score("22", "11", [[p["L_LOH"]]]) == 1              # one interval loss
assert o.score("21", "10", [[p["L_LOH"]]]) == 1              # one interval loss including zero
assert o.score("2X2", "1X1", [[p["L_LOH"]]]) == 2            # X breaks intervals
A = [p["L_LOH"], p["G"]]
assert math.isinf(o.score("0", "1", [A]))                    # 0 -> 1: gains cannot resurrect a lost allele
assert o.score("1", "0", [A]) == 1 and o.score("1", "2", [A]) == 1
assert o.pair("0", "1", [A], [A]) == 1                       # shared ancestor, zero WGDs: 1 loses one copy toward 0
# the full edge model P = A + ADA + ADADA: at most two WGDs
assert o.directed("1X1", "8X8") == 4                         # 2 gains + 2 WGDs
assert math.isinf(o.directed("0", "1")) and o.kernel("0", "1") == 1
assert o.kernel("0X0", "7X7") == 8 and o.kernel("1X1", "8X8") == o.kernel("8X8", "1X1")
assert o.optimum({"1": 0}, {"2": 0}) == {"1": 1, "2": 1}       # both tied ancestors survive
assert o.optimum({"11": 0}, {"22": 0}) == {"11": 1, "22": 1}   # NOT {"12","21"}: correlations kept
# P equals the chain A (W1 A)^2 with MEDICC's optional capped WGD W1 (Step 1)
W1 = medicc.factory.create_1step_WGD_fst(o.sym, "X", wgd_cost=1, minimize=True, wgd_x2=True, wgd_cap=True)
optional = [A + [W1, *A] * 2]
for x in ("12", "0X3", "28"):                                 # every child word, with its cost
    assert _words(o.children(o.language({x: 0})), o.sym) == _words(o.children(o.language({x: 0}), optional), o.sym)
print("PASS step 2")
