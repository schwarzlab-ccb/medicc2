from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import fstlib
from medicc import io
from medicc.fstengine import ModelConfig
assert io.read_fst() == ModelConfig()
assert io.read_fst(no_wgd=True, total_copy_numbers=True) == ModelConfig(True, True)
assert isinstance(io.read_fst(user_fst=str(ROOT / "medicc/objects/wgd_asymm.fst")), fstlib.Fst)
print("PASS step 15")