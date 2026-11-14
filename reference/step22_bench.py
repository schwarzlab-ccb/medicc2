"""Phase 5 benchmark: time the fstengine on meso5 with 1 core, and check the answers never change.

REFERENCE for migration Step 22 (not part of REF; untracked). Goes to tests/migration/bench.py.
Written for this plan: REF's tools/benchmark_native.py (the closest match) times the finished
engine for PGO builds on gundem data, with NNI and events; none of that exists yet here.

Three timings (best of --repeat runs):
  pairs   the full pairwise distance matrix (medicc.fstengine.pairwise_distances, 1 core)
  slowest the slowest single pair, timed on its own with a fresh Model
  tree    one reconstruction: medicc.main with the fixed tree medicc/test_data/meso5.new
          (meso5 only; for another --input medicc.main builds the NJ tree, distances included)
The Phase 4 engine needs minutes for ONE full meso5 pair, so the default input is meso5 cut to
its first 4 chromosomes. After Step 23 run all of meso5 (--chromosomes 0); after Step 26 meso5
takes about 0.01 s, so time a harder input too, e.g. gundem PTX010.

The first run of an input saves its answers (distances, tree length, ancestral profiles) to
tests/migration/bench_answers_<input>_<size>.json; later runs must give the same distances and
tree length. The ancestral profiles may change to another optimal tie: that only prints a note.

    python -B tests/migration/bench.py                    # meso5, first 4 chromosomes
    python -B tests/migration/bench.py --chromosomes 0    # all of meso5
    python -B tests/migration/bench.py --chromosomes 0 --input examples/gundem_et_al_2015/PTX010_input_df.tsv
"""
import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from Bio import Phylo
import medicc
from medicc.fstengine import ModelConfig, pairwise_distances
from medicc.fstengine.model import Model


MESO5 = ROOT / 'medicc/test_data/meso5.tsv'


def load(path, chromosomes):
    """The input, cut to its first `chromosomes` chromosomes (0 = all)."""
    df = medicc.io.read_and_parse_input_data(path)
    if chromosomes:
        keep = df.index.get_level_values('chrom').unique()[:chromosomes]
        df = df.loc[df.index.get_level_values('chrom').isin(keep)]
    return df


def time_pairs(profiles):
    """The distance matrix, as MEDICC computes it (scheduler, 1 core).
    Returns (seconds, distances as {"x:y": d})."""
    start = time.perf_counter()
    matrix = pairwise_distances(profiles, ModelConfig(), n_cores=1)
    seconds = time.perf_counter() - start
    names = list(profiles)
    return seconds, {f'{x}:{y}': float(matrix.loc[x, y]) for i, x in enumerate(names) for y in names[i + 1:]}


def slowest_pair(profiles):
    """Every pair timed on its own, one fresh Model each (no cache help). Returns (seconds, pair)."""
    names = list(profiles)
    slowest = (0.0, None)
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            with Model() as model:
                start = time.perf_counter()
                model.score(profiles[x], profiles[y])
                slowest = max(slowest, (time.perf_counter() - start, f'{x}:{y}'))
    return slowest


def time_tree(df, fixed_tree):
    """One reconstruction. Returns (seconds, tree length, hash of the ancestral profiles)."""
    tree = Phylo.read(ROOT / 'medicc/test_data/meso5.new', 'newick') if fixed_tree else None
    start = time.perf_counter()
    *_, final_tree, output_df, _ = medicc.main(df, ModelConfig(), input_tree=tree, n_cores=1)
    seconds = time.perf_counter() - start
    return seconds, final_tree.total_branch_length(), hashlib.sha256(output_df.to_csv().encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--chromosomes', type=int, default=4, help='first N chromosomes; 0 = all')
    parser.add_argument('--input', type=Path, default=MESO5, help='MEDICC input .tsv (default: meso5)')
    parser.add_argument('--repeat', type=int, default=3, help='runs per timing; the best is printed')
    args = parser.parse_args()
    logging.disable(logging.CRITICAL)
    warnings.filterwarnings('ignore', category=FutureWarning)

    path = args.input if args.input.is_absolute() else ROOT / args.input
    df = load(path, args.chromosomes)
    profiles = medicc.create_cn_str_dict(df)
    name, size = path.name.split('.')[0].replace('_input_df', ''), args.chromosomes or 'all'
    print(f'{name}, chromosomes: {size}, {len(profiles)} samples, {len(next(iter(profiles.values())))} columns')

    runs = [time_pairs(profiles) for _ in range(args.repeat)]
    pair_seconds, distances = min(r[0] for r in runs), runs[0][1]
    slow_seconds, slow_pair = min(slowest_pair(profiles) for _ in range(args.repeat))
    print(f'pairs   {pair_seconds:8.3f} s  ({len(distances)} pairs)')
    print(f'slowest {slow_seconds:8.3f} s  ({slow_pair}, distance {distances[slow_pair]:g})')
    runs = [time_tree(df, path == MESO5) for _ in range(args.repeat)]
    tree_seconds, length, profile_hash = min(r[0] for r in runs), runs[0][1], runs[0][2]
    print(f'tree    {tree_seconds:8.3f} s  (tree length {length:g})')

    answers = {'distances': distances, 'tree_length': length, 'profiles_sha256': profile_hash}
    saved = Path(__file__).with_name(f'bench_answers_{name}_{size}.json')
    if saved.exists():
        before = json.loads(saved.read_text())
        if (before['distances'], before['tree_length']) != (distances, length):
            sys.exit(f'FAIL: distances or tree length differ from {saved.name}')
        print(f'same distances and tree length as {saved.name}')
        if before['profiles_sha256'] != profile_hash:
            # Same tree length, other profiles: a different choice among equally good ancestors
            # (e.g. Step 27 changes the order of tied states). Not an error; say so in the commit.
            print('note: other ancestral profiles of the same cost (another tie was chosen)')
    else:
        saved.write_text(json.dumps(answers, indent=2) + '\n')
        print(f'saved answers to {saved.name}')


if __name__ == '__main__':
    main()
