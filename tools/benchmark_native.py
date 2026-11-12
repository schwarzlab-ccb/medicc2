#!/usr/bin/env python3
"""Train PGO or benchmark the currently built engine on held-out bundled data.

Run from a source checkout. Training and benchmark datasets are disjoint.
Calibration fixes iteration counts for subsequent standard/PGO comparisons.
Imports, input reading, warmup and output hashing are outside measured regions;
model creation and all computation in each operation are included.
"""
import argparse
import gc
import hashlib
import json
import logging
import math
from pathlib import Path
import sys
import time
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import medicc
from medicc.native import ModelConfig, pairwise_distances


def digest(result):
    """Hash distances, topology/lengths, profiles and events, not just a total score."""
    if hasattr(result, 'to_csv'):
        payload = result.to_csv()
    else:
        labels, distances, _, tree, profiles, events = result
        def clade(node):
            return [node.name, node.branch_length, [clade(c) for c in node.clades]]
        payload = json.dumps([list(labels), clade(tree.root)], sort_keys=True)
        for frame in (distances, profiles, events):
            payload += '\n' + (frame.to_csv() if frame is not None else 'None')
    return hashlib.sha256(payload.encode()).hexdigest()


def datasets():
    return {patient: medicc.io.read_and_parse_input_data(
        ROOT / f'examples/gundem_et_al_2015/{patient}_input_df.tsv')
        for patient in ('PTX004', 'PTX007', 'PTX010')}


def workloads():
    frames = datasets()
    tasks = {}
    for patient, frame in frames.items():
        profiles = medicc.create_cn_str_dict(frame)
        tasks[f'pairwise_{patient}'] = (
            lambda profiles=profiles: pairwise_distances(profiles, ModelConfig(), n_cores=1))
    for patient in ('PTX004', 'PTX007'):
        frame = frames[patient]
        tasks[f'pipeline_events_{patient}'] = (
            lambda frame=frame: medicc.main(frame, ModelConfig(), n_cores=1,
                                           tree_method='nj', reconstruct_events=True))
    frame = medicc.io.read_and_parse_input_data(ROOT / 'medicc/test_data/meso5.tsv')
    tasks['pipeline_nni_meso5'] = lambda: medicc.main(
        frame, ModelConfig(), n_cores=1, tree_method='nj', nni=True,
        reconstruct_events=True)
    return tasks


def train():
    """Exercise pair scoring, joins, resolution, histories, NNI and both data modes."""
    for relative in ('examples/simple_example/simple_example.tsv',
                     'examples/testing_example/testing_example.tsv',
                     'examples/gundem_et_al_2015/PTX005_input_df.tsv'):
        frame = medicc.io.read_and_parse_input_data(ROOT / relative)
        for no_wgd in (False, True):
            medicc.main(frame, ModelConfig(no_wgd=no_wgd), n_cores=2,
                        tree_method='nj', reconstruct_events=True)
        print(f'Trained {relative}', flush=True)
    frame = medicc.io.read_and_parse_input_data(ROOT / 'examples/simple_example/simple_example.tsv')
    medicc.main(frame, ModelConfig(), n_cores=1, tree_method='nj', nni=True,
                reconstruct_events=True)
    frame = medicc.io.read_and_parse_input_data(
        ROOT / 'examples/simple_example/simple_example.tsv',
        allele_columns=['cn_a'], total_copy_numbers=True)
    medicc.main(frame, ModelConfig(total_copy_numbers=True), n_cores=1,
                allele_columns=['cn_a'], total_cn=True, tree_method='nj', reconstruct_events=True)
    print('Training complete; profiles are flushed when this process exits.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--train', action='store_true')
    parser.add_argument('--calibrate', type=Path, help='Write fixed iteration counts as JSON')
    parser.add_argument('--seconds', type=float, default=1.0, help='Target seconds per calibrated workload')
    parser.add_argument('--iterations', type=Path, help='Read calibrated iteration counts')
    parser.add_argument('--output', type=Path, help='Write timings and output hashes as JSON')
    parser.add_argument('--workload', action='append', help='Select named workloads; default: all six')
    args = parser.parse_args()
    if args.seconds <= 0:
        parser.error('--seconds must be positive')
    logging.disable(logging.CRITICAL)
    warnings.filterwarnings('ignore', category=FutureWarning)
    if args.train:
        train()
        return
    tasks = workloads()
    selected = args.workload or list(tasks)
    if set(selected) - tasks.keys():
        parser.error('Unknown workload: ' + ', '.join(set(selected) - tasks.keys()))
    counts = json.loads(args.iterations.read_text()) if args.iterations else {}
    results = {}
    for name in selected:
        operation = tasks[name]
        # A fresh model is created for every call, so warmup does not retain its caches.
        expected = digest(operation())
        start = time.perf_counter()
        sample = operation()
        sample_time = time.perf_counter() - start
        if digest(sample) != expected:
            raise RuntimeError(f'Non-reproducible warmup output in {name}')
        count = int(counts.get(name, max(1, min(1000, math.ceil(args.seconds / sample_time)))))
        if count < 1:
            parser.error(f'Iteration count for {name} must be positive')
        counts[name] = count
        timings = []
        for _ in range(count):
            gc.collect()
            start = time.perf_counter()
            result = operation()
            timings.append(time.perf_counter() - start)
            if digest(result) != expected:
                raise RuntimeError(f'Non-reproducible output in {name}')
        results[name] = dict(iterations=count, seconds_per_operation=float(np.mean(timings)),
                             seconds=timings, output_sha256=expected)
        print(f'{name}: {results[name]["seconds_per_operation"]:.6f}s ({count} iterations)', flush=True)
    if args.calibrate:
        args.calibrate.write_text(json.dumps(counts, indent=2) + '\n')
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
