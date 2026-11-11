"""Tests of the native OpenFST model (medicc.native), the default model of MEDICC2."""
import ctypes
import io
import os
import pathlib
import pickle
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

import medicc
from medicc.native import Model, ModelConfig
from medicc.native.model import PARTS_FILE, Message

ROOT = pathlib.Path(__file__).parent.parent.absolute()
DATA = pathlib.Path(__file__).parent / 'test_data'
INF = 100000000


@pytest.fixture(scope='module')
def two_bin_oracle():
    """All 81 two-bin profiles with directed and undirected distances from an independent
    enumeration of the literal L_LOH, G and capped doubling primitives."""
    fixture = np.load(DATA / 'two_bin_oracle.npz')
    return fixture['words'].tolist(), fixture['directed'], fixture['undirected']


@pytest.fixture(scope='module')
def a_only_oracle():
    """Directed and undirected two-bin distances of A = L_LOH o G, computed by a literal
    DP over the primitives in objects/ancestral_parts.txt (no native composition)."""
    words = [f"{a}{b}" for a in range(9) for b in range(9)]
    tokens = iter(pathlib.Path(PARTS_FILE).read_text().split())
    relations = {}
    for name in tokens:
        n, start, count = (int(next(tokens)) for _ in range(3))
        finals = np.array([int(next(tokens)) for _ in range(n)], dtype=float)
        finals[finals >= INF] = np.inf
        arcs = [tuple(int(next(tokens)) for _ in range(5)) for _ in range(count)]
        if name not in ('L_LOH', 'G'):
            continue
        by_labels = {}
        for source, a, b, cost, target in arcs:
            by_labels.setdefault((a, b), []).append((source, cost, target))
        matrix = np.full((len(words), len(words)), np.inf)
        for i, left in enumerate(words):
            for j, right in enumerate(words):
                costs = np.full(n, np.inf)
                costs[start] = 0
                for a, b in zip(left, right):
                    step = np.full(n, np.inf)
                    for source, cost, target in by_labels.get((int(a), int(b)), []):
                        step[target] = min(step[target], costs[source] + cost)
                    costs = step
                matrix[i, j] = min(costs + finals)
        relations[name] = matrix
    directed = np.min(relations['L_LOH'][:, :, None] + relations['G'][None, :, :], axis=1)
    undirected = np.min(directed[:, :, None] + directed[:, None, :], axis=0)
    return words, directed, undirected


def weighted_message(model, words):
    """Build a message that accepts each word at the given cost."""
    rows = [[INF] + [-1, 0] * 10]
    for word, cost in words.items():
        state = 0
        for char in word:
            column = 1 + 2 * (9 if char == 'X' else int(char))
            if rows[state][column] < 0:
                rows[state][column] = len(rows)
                rows.append([INF] + [-1, 0] * 10)
            state = rows[state][column]
        rows[state][0] = cost
    values = (ctypes.c_int * (21 * len(rows)))(*(v for row in rows for v in row))
    return Message(model, model.lib.am_import(len(next(iter(words))), len(rows), 0, 0, values))


def language(message):
    """Enumerate all words of a (small) message with their costs."""
    message.check_open()
    lib, handle = message.model.lib, message.handle
    result = {}
    stack = [(lib.am_start(handle), '', lib.am_offset(handle))]
    while stack:
        state, prefix, cost = stack.pop()
        if state < 0:
            continue
        final = lib.am_final(handle, state)
        if final < INF:
            result[prefix] = cost + final
        for c in range(10):
            weight = ctypes.c_int()
            target = lib.am_arc(handle, state, c, ctypes.byref(weight))
            if target >= 0:
                stack.append((target, prefix + ('X' if c == 9 else str(c)), cost + weight.value))
    return result


def branch_frame(parent, child, column='cn'):
    """Input DataFrame with the samples 'diploid' (parent) and 'sample' (child)."""
    rows = []
    for name, text in (('diploid', parent), ('sample', child)):
        for c, word in enumerate(text.split('X'), 1):
            for i, value in enumerate(word):
                rows.append((name, 'chr' + str(c), i * 10, i * 10 + 10, value))
    df = pd.DataFrame(rows, columns=['sample_id', 'chrom', 'start', 'end', column])
    df['chrom'] = medicc.tools.format_chromosomes(df['chrom'])
    return df.set_index(['sample_id', 'chrom', 'start', 'end']).sort_index()


def run_medicc2(*args):
    return subprocess.run([sys.executable, 'medicc2', *map(str, args)], capture_output=True, text=True,
                          cwd=ROOT, timeout=600, env={**os.environ, 'MPLBACKEND': 'Agg'})


def test_two_bin_distances_match_oracle(two_bin_oracle):
    words, directed, undirected = two_bin_oracle
    with Model() as model:
        for i, x in enumerate(words):
            for j, y in enumerate(words):
                assert model.score(x, y) == undirected[i, j], (x, y)
                assert model.score(x, y, directed=True) == directed[i, j], (x, y)
        model.clear()
        assert model.score('01', '88') == undirected[1, 80]


def test_joins_keep_all_minimum_cost_ancestors(two_bin_oracle):
    words, directed, _ = two_bin_oracle
    rng = np.random.default_rng(1130)
    with Model() as model:
        for _ in range(100):
            left = {words[i]: int(rng.integers(12)) for i in rng.choice(81, 5, False)}
            right = {words[i]: int(rng.integers(12)) for i in rng.choice(81, 5, False)}
            a = np.array([left.get(w, np.inf) for w in words])
            b = np.array([right.get(w, np.inf) for w in words])
            costs = np.min(directed + a[None, :], axis=1) + np.min(directed + b[None, :], axis=1)
            expected = {w: int(c) for w, c in zip(words, costs) if c == costs.min()}
            x, y = weighted_message(model, left), weighted_message(model, right)
            z = model.join(x, y)
            assert language(z) == expected
            for message in (x, y, z):
                message.close()


def test_resolve_with_correlations_and_large_weights():
    languages = [{'2X4': 0, '4X8': 1}, {'2XX4': 3, '4XX4': 4}, {'X2X4': 0, 'X4X4': 2},
                 {'2X': 0, '4X': 2}, {'02': 2**24 + 3, '22': 2**24 + 4}]
    with Model() as model:
        for words in languages:
            message = weighted_message(model, words)
            for parent in words:
                child, cost = model.resolve(parent, message)
                assert cost == min(model.score(parent, w, True) + v for w, v in words.items())
                assert model.score(parent, child, True) + words[child] == cost
            message.close()


def test_validation_and_handle_lifetimes():
    with pytest.raises(ValueError):
        ModelConfig(no_wgd='yes')
    with Model() as model:
        with pytest.raises(ValueError):
            model.score('1X', '11')
        with pytest.raises(ValueError):
            model.score('9', '1')
        message = model.string('2')
    # messages own their FST independently of the model handle
    assert language(message) == {'2': 0}
    message.close()
    message.close()
    with pytest.raises(ValueError):
        language(message)


def test_total_copy_numbers_double_instead_of_adding_two(two_bin_oracle):
    words, directed, _ = two_bin_oracle
    config = medicc.io.read_fst(total_copy_numbers=True)
    assert config.total_copy_numbers
    with Model(config) as total, Model() as allele:
        for x, y, expected in [('3X3X3X3', '6X6X6X6', 1), ('2X2X2X2', '8X8X8X8', 2),
                               ('5X5X5X5', '8X8X8X8', 1)]:
            assert total.score(x, y, True) == expected
            assert total.score(x, y) == allele.score(x, y)
        for i in (0, 1, 9, 13, 22, 40, 80):
            for j in (0, 5, 17, 53, 80):
                assert total.score(words[i], words[j], True) == directed[i, j]


def test_no_wgd_distances_match_a_only_oracle(a_only_oracle):
    words, directed, undirected = a_only_oracle
    with Model(ModelConfig(no_wgd=True)) as model:
        for i, x in enumerate(words):
            for j, y in enumerate(words):
                assert model.score(x, y, True) == directed[i, j], (x, y)
                assert model.score(x, y) == undirected[i, j], (x, y)
        assert model.score('1X1X1X1', '4X4X4X4', True) == 12
        assert model.score('2X2X2X2', '8X8X8X8', True) == 24


def test_no_wgd_joins_and_resolution_use_a_only(a_only_oracle):
    words, directed, _ = a_only_oracle
    rng = np.random.default_rng(1111)
    with Model(ModelConfig(no_wgd=True)) as model:
        for _ in range(30):
            left, right = [{words[i]: int(rng.integers(8)) for i in rng.choice(81, 5, False)}
                           for _ in range(2)]
            a, b = [np.array([d.get(w, np.inf) for w in words]) for d in (left, right)]
            costs = np.min(directed + a, axis=1) + np.min(directed + b, axis=1)
            expected = {w: c for w, c in zip(words, costs) if c == min(costs)}
            x, y = weighted_message(model, left), weighted_message(model, right)
            z = model.join(x, y)
            assert language(z) == expected
            i = int(rng.choice([10, 20, 40, 80]))
            child, cost = model.resolve(words[i], x)
            assert cost == min(directed[i] + a)
            assert cost == directed[i, words.index(child)] + left[child]
            for message in (x, y, z):
                message.close()


def test_event_histories_replay_and_count_equals_branch_length(two_bin_oracle):
    words, directed, _ = two_bin_oracle
    rng = np.random.default_rng(733)
    pairs = [('213', '000'), ('101', '202'), ('2X2X2X2', '8X8X8X8'), ('0X11', '0X11')]
    pairs += [(''.join(map(str, rng.integers(1, 9, 3))), ''.join(map(str, rng.integers(0, 9, 3))))
              for _ in range(40)]
    pairs += [(words[i], words[j]) for i, j in rng.integers(81, size=(35, 2))
              if np.isfinite(directed[i, j])]
    for config in (ModelConfig(), ModelConfig(no_wgd=True)):
        with Model(config) as model:
            for parent, child in pairs:
                score = model.score(parent, child, True)
                if not np.isfinite(score):
                    continue
                tree = Phylo.read(io.StringIO(f"(diploid:0,sample:{score});"), 'newick')
                frame = branch_frame(parent, child)
                profiles, events = medicc.event_reconstruction.calculate_all_cn_events(
                    tree, frame, ['cn'], model_config=config)
                assert len(events) == score, (parent, child)
                # replay the exported intervals independently of the native stage strings
                current = frame.loc['diploid', 'cn'].astype(int).copy()
                for event in events.reset_index().sort_values('event_order').itertuples():
                    if event.type == 'wgd':
                        current = (current * 2).clip(upper=8)
                        continue
                    idx = [k for k in current.index
                           if k[0] == event.chrom and k[1] >= event.start and k[2] <= event.end]
                    values = current.loc[idx]
                    current.loc[idx] = ((values - 1).clip(lower=0) if event.type == 'loss'
                                        else values.where(values == 0, values + 1))
                assert current.tolist() == frame.loc['sample', 'cn'].astype(int).tolist()
                assert profiles.loc['sample', 'cn'].tolist() == current.tolist()


def test_main_with_fixed_tree_and_cores():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    fst = medicc.io.read_fst()
    for n_cores in (1, 2):
        tree = medicc.io.import_tree(DATA / 'meso5.new')
        labels, distances, _, final_tree, output_df, events_df = medicc.main(
            input_df, fst, input_tree=tree, n_cores=n_cores, reconstruct_events=True)
        assert final_tree.total_branch_length() == 67
        assert (distances.values == 0).all()
        assert len(events_df) == 67
        for node in final_tree.find_clades():
            if node.name not in (None, 'diploid'):
                assert (events_df.index.get_level_values('sample_id') == node.name).sum() == \
                    node.branch_length
        for sample in labels:
            pd.testing.assert_frame_equal(output_df.loc[sample, ['cn_a', 'cn_b']].astype(str),
                                          input_df.loc[sample], check_names=False)


@pytest.mark.fastme
def test_main_topology_only_and_single_tumour():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    fst = medicc.io.read_fst()
    labels, distances, _, tree, _, events_df = medicc.main(input_df, fst, ancestral_reconstruction=False,
                                                           n_cores=2)
    assert distances.loc['diploid', 'leaf_0'] == 19
    assert distances.loc['leaf_2', 'leaf_4'] == 41
    assert events_df is None
    assert {c.name for c in tree.get_terminals()} == set(labels)

    two = input_df.loc[['diploid', 'leaf_0']]
    _, distances, _, tree, _, _ = medicc.main(two, fst)
    assert distances.loc['diploid', 'leaf_0'] == 19
    assert tree.total_branch_length() == 19


@pytest.mark.parametrize('total_cn', [False, True])
def test_main_no_wgd(total_cn):
    normal, child = ('2', '8') if total_cn else ('1', '4')
    frame = branch_frame('X'.join([normal] * 4), 'X'.join([child] * 4))
    for fst, no_wgd in [(medicc.io.read_fst(total_copy_numbers=total_cn, no_wgd=True), False),
                        (medicc.io.read_fst(), True)]:
        result = medicc.main(frame, fst, no_wgd=no_wgd, total_cn=total_cn, allele_columns=['cn'],
                             reconstruct_events=True, n_cores=2)
        assert result[3].total_branch_length() == (24 if total_cn else 12)
        assert len(result[5]) == result[3].total_branch_length()
        assert set(result[5]['type']) == {'gain'}
    with_wgd = medicc.main(frame, medicc.io.read_fst(), allele_columns=['cn'],
                           reconstruct_events=True)
    assert 'wgd' in set(with_wgd[5]['type'])
    assert with_wgd[3].total_branch_length() == len(with_wgd[5])


def test_main_rejects_unsupported_native_options():
    frame = branch_frame('11', '22')
    fst = medicc.io.read_fst()
    with pytest.raises(medicc.core.MEDICCError, match='work in progress'):
        medicc.main(frame, fst, prune_weight=1)
    with pytest.raises(medicc.core.MEDICCError, match='one copy-number column'):
        medicc.main(frame.assign(cn_b=frame['cn']), fst, total_cn=True)
    with pytest.raises(TypeError):
        medicc.io.read_fst(n_wgd=1)


def test_read_fst_returns_native_model():
    assert medicc.io.read_fst() == ModelConfig()
    assert medicc.io.read_fst(no_wgd=True, total_copy_numbers=True) == \
        ModelConfig(total_copy_numbers=True, no_wgd=True)


@pytest.mark.fastme
def test_user_fst_reproduces_master():
    """A custom FST is loaded with fstlib (no_wgd is ignored) and master's bundled FSTs give
    master's results (tree length 608 / 683); events are rejected for a custom FST."""
    input_df = medicc.io.read_and_parse_input_data(ROOT / 'examples/testing_example/testing_example.tsv')
    objects = ROOT / 'medicc' / 'objects'
    for name, expected in (('wgd_asymm.fst', 608), ('no_wgd_asymm.fst', 683)):
        fst = medicc.io.read_fst(user_fst=str(objects / name), no_wgd=True)
        reference = medicc.io.read_fst(user_fst=str(objects / name))
        result = medicc.main(input_df, fst)
        assert result[3].total_branch_length() == expected
        assert result[1].equals(medicc.main(input_df, reference)[1])
    with pytest.raises(medicc.core.MEDICCError, match='event labels'):
        medicc.main(input_df, reference, reconstruct_events=True)


@pytest.mark.fastme
def test_native_no_wgd_equals_master_no_wgd_tree_length():
    input_df = medicc.io.read_and_parse_input_data(ROOT / 'examples/testing_example/testing_example.tsv')
    result = medicc.main(input_df, medicc.io.read_fst(no_wgd=True), n_cores=2)
    assert result[3].total_branch_length() == 683


@pytest.mark.fastme
def test_bootstrap_with_native_model():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    fst = medicc.io.read_fst()
    tree = medicc.main(input_df, fst)[3]
    trees_df, support_tree = medicc.bootstrap.run_bootstrap(input_df, original_tree=tree, N_bootstrap=5,
                                                            seed=1, fst=fst, show_progress=False,
                                                            n_cores=2)
    assert trees_df['count'].sum() == 5
    assert support_tree is not None
    pickle.dumps(fst)


@pytest.mark.parametrize('n_bootstrap', [0, 1, 4])
def test_bootstrap_rejects_fewer_than_five_runs(n_bootstrap):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    with pytest.raises(medicc.core.MEDICCError, match='at least 5'):
        medicc.bootstrap.run_bootstrap(input_df, N_bootstrap=n_bootstrap, show_progress=False)


def test_overlap_events_keeps_every_wgd():
    events = pd.DataFrame(
        [['internal_1', 'chr0', 0, 100, 'both', 'wgd', 0],
         ['leaf_1', 'chr0', 0, 100, 'both', 'wgd', 0],
         ['leaf_1', 'chr0', 0, 100, 'both', 'wgd', 0],
         ['leaf_1', 'chr1', 0, 1000000, 'cn_a', 'gain', 2]],
        columns=medicc.event_reconstruction.EVENT_COLUMNS).set_index(['sample_id', 'chrom', 'start', 'end'])
    overlaps = medicc.event_reconstruction.overlap_events(events_df=events, regions_bed=None)
    wgd = overlaps.loc[overlaps['event'] == 'WGD']
    assert sorted(wgd.index) == ['internal_1', 'leaf_1', 'leaf_1']
    assert (wgd['Chromosome'] == 'all').all()
    assert (wgd[['Event-Start', 'Event-End']] == 0).all().all()
    assert (wgd['final_name'] == 'WGD').all()


@pytest.mark.fastme
def test_cli_events_beds_plots_and_bootstrap(tmp_path):
    frame = branch_frame('2X2X2X2', '8X8X8X8')
    extra = frame.loc[['sample']].rename(index={'sample': 'other'}).replace('8', '6')
    pd.concat([frame, extra]).to_csv(tmp_path / 'input.tsv', sep='\t')
    (tmp_path / 'regions.bed').write_text('chr1\t0\t10\tREGION\n')
    outputs = {}
    for name, flags in (('default', []), ('wgd_x2', ['--wgd-x2'])):
        out = tmp_path / name
        run = run_medicc2(tmp_path / 'input.tsv', out, '--total-copy-numbers', '-a', 'cn',
                          '--events', '--regions-bed', tmp_path / 'regions.bed', '--chromosomes-bed',
                          'default', '--plot', 'both', '--prefix', 'case', '--bootstrap-nr', '5',
                          '-j', '2', *flags)
        assert run.returncode == 0, run.stderr
        assert ('--wgd-x2 has no effect' in run.stderr) == (name == 'wgd_x2')
        events = pd.read_csv(out / 'case_copynumber_events_df.tsv', sep='\t')
        summary = pd.read_csv(out / 'case_summary.tsv', sep='\t', header=None, index_col=0)[1]
        assert len(events) == float(summary['tree_length'])
        assert 'wgd' in set(events['type'])
        profiles = pd.read_csv(out / 'case_final_cn_profiles.tsv', sep='\t')
        assert set(profiles.loc[profiles['sample_id'] == 'diploid', 'cn']) == {2}
        for suffix in ('_events_overlap.tsv', '_cn_profiles.pdf', '_cn_profiles_heatmap.pdf',
                       '_final_tree.png', '_support_tree.pdf', '_support_tree.new',
                       '_bootstrap_trees_df.pickle', '_pairwise_distances.tsv', '_branch_lengths.tsv'):
            assert (out / ('case' + suffix)).is_file(), suffix
        outputs[name] = (out / 'case_copynumber_events_df.tsv').read_bytes()
    assert outputs['default'] == outputs['wgd_x2']


@pytest.mark.fastme
def test_cli_tree_topology_only_and_no_wgd(tmp_path):
    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'tree', '--tree', DATA / 'meso5.new',
                      '--plot', 'none', '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    assert Phylo.read(tmp_path / 'tree/case_final_tree.new', 'newick').total_branch_length() == 67

    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'topology', '--topology-only',
                      '--plot', 'heatmap', '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    assert not (tmp_path / 'topology/case_final_cn_profiles.tsv').exists()
    distances = pd.read_csv(tmp_path / 'topology/case_pairwise_distances.tsv', sep='\t', index_col=0)
    assert distances.loc['leaf_2', 'leaf_4'] == 41

    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'nowgd', '--no-wgd', '--events',
                      '--plot', 'none', '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    events = pd.read_csv(tmp_path / 'nowgd/case_copynumber_events_df.tsv', sep='\t')
    assert 'wgd' not in set(events['type'])


@pytest.mark.parametrize('flags, message', [
    (['--prune-weight', '1'], 'work in progress'),
    (['--fst', 'medicc/objects/wgd_asymm.fst', '--events'], 'event labels'),
    (['--bootstrap-nr', '4'], 'at least 5'),
    (['--fst-chr-separator', 'Y'], 'unrecognized arguments'),
    (['--input-chr-separator', 'Y'], 'unrecognized arguments'),
    (['--native'], 'unrecognized arguments')])
def test_cli_rejects_unsupported_options(tmp_path, flags, message):
    out = tmp_path / 'invalid'
    run = run_medicc2(DATA / 'meso5.tsv', out, *flags)
    assert run.returncode == 2
    assert message in run.stderr
    assert not out.exists()


@pytest.mark.fastme
def test_cli_fst_reproduces_master_tree_length(tmp_path):
    out = tmp_path / 'fst'
    run = run_medicc2('examples/testing_example/testing_example.tsv', out, '--fst',
                      'medicc/objects/wgd_asymm.fst', '--plot', 'none', '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    summary = pd.read_csv(out / 'case_summary.tsv', sep='\t', header=None, index_col=0)[1]
    assert float(summary['tree_length']) == 608


def test_write_branch_lengths_leaves_out_the_given_normal(tmp_path):
    tree = Phylo.read(io.StringIO("((a:2,b:3)internal1:1,my_normal:0)root:0;"), 'newick')
    out = tmp_path / 'branch_lengths.tsv'
    medicc.io.write_branch_lengths(tree, out, normal_name='my_normal')
    names = [line.split('\t')[0] for line in out.read_text().splitlines()]
    assert names == ['root', 'internal1', 'a', 'b']
