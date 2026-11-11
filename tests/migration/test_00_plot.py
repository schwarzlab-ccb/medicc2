from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import warnings
import matplotlib
matplotlib.use("Agg")
import medicc

frame = medicc.io.read_and_parse_input_data(ROOT / "examples/simple_example/simple_example.tsv")
with warnings.catch_warnings():
    warnings.simplefilter("error")
    medicc.plot.plot_cn_profiles(frame)
    medicc.plot.plot_cn_heatmap(frame)
print("PASS step 0.3")
