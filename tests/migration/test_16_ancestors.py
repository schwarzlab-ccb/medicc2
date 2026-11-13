from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
from io import StringIO
from Bio import Phylo
from medicc.ancestors import reconstruct_ancestors
from medicc.fstengine import ModelConfig
tree = Phylo.read(StringIO("(diploid,a,b)r;"), "newick")
result = reconstruct_ancestors(tree, {"diploid": "11", "a": "11", "b": "22"}, ModelConfig(), "diploid", n_cores=2)
assert result["r"] == "11" and result["b"] == "22"
print("PASS step 16")