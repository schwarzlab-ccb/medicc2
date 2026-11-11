from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import math
import fstlib
import medicc

sym = medicc.create_symbol_table(8, "X")
fsa = lambda s: fstlib.factory.from_string(s, isymbols=sym, osymbols=sym)

def score(W, a, b):
    return float(fstlib.score(W, fsa(a), fsa(b)))

capped = medicc.factory.create_1step_WGD_fst(sym, "X", wgd_cost=1, minimize=True, wgd_x2=True, wgd_cap=True)
plain = medicc.factory.create_1step_WGD_fst(sym, "X", wgd_cost=1, minimize=True, wgd_x2=True)
assert score(capped, "5", "8") == 1          # saturation: 5 -> min(10, 8)
assert math.isinf(score(plain, "5", "8"))    # master behaviour unchanged
assert score(capped, "4", "8") == 1
assert score(capped, "2X3", "4X6") == 1      # one WGD for the whole genome
assert score(capped, "3", "3") == 0          # the WGD is optional in this 1-step FST
try:
    medicc.factory.create_1step_WGD_fst(sym, "X", wgd_cap=True)
except ValueError:
    pass
else:
    raise AssertionError("wgd_cap without wgd_x2 accepted")
print("PASS step 1")
