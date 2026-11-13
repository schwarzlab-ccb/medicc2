from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import medicc
frame = medicc.io.read_and_parse_input_data(ROOT / "medicc/test_data/meso5.tsv")
tree = medicc.io.import_tree(ROOT / "medicc/test_data/meso5.new")
out = medicc.main(frame, medicc.io.read_fst(), input_tree=tree, n_cores=2, reconstruct_events=False)
assert out[3].total_branch_length() == 67
print("PASS step 17")