"""Tests of the NNI topology refinement (medicc.nni, --nni-mode)."""
import io
import os
import pathlib
import subprocess
import sys

import pandas as pd
import pytest
from Bio import Phylo

import medicc
from medicc.nni import nni_neighbours, nni_search, topology_key

ROOT = pathlib.Path(__file__).parent.parent.absolute()
DATA = pathlib.Path(__file__).parent / 'test_data'

# meso5 with shuffled leaves: tree length 98, NNI recovers the tree of meso5.new (67)
BAD_TREE = "(diploid,((leaf_1,(leaf_3,leaf_0)internal_4)internal_2,(leaf_4,leaf_2)internal_3)internal_1);"


def import_tree(newick):
    """Tree as returned by medicc.io.import_tree (input_tree of medicc.main)."""
    return Phylo.read(io.StringIO(newick), 'newick')


def read_tree(newick):
    """Tree rooted as in medicc.main for ancestor reconstruction (normal is a child of the root)."""
    tree = import_tree(newick)
    tree.root_with_outgroup([c for c in tree.root.clades if c.name != 'diploid'][0].name)
    return tree


@pytest.fixture(scope='module')
def meso5():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    cn_str_dict = medicc.create_cn_str_dict(input_df)
    fst = medicc.io.read_fst()

    def evaluate(tree, n_cores=1):
        return medicc.reconstructed_tree_length(tree, cn_str_dict, fst, 'diploid', n_cores=n_cores)

    return input_df, evaluate


def run_medicc2(*args):
    return subprocess.run([sys.executable, 'medicc2', *map(str, args)], capture_output=True, text=True,
                          cwd=ROOT, timeout=600, env={**os.environ, 'MPLBACKEND': 'Agg'})


def test_moves_two_per_internal_edge_distinct_and_names_kept():
    # internal edges below the root with an internal child: (i0,i1), (i0,i2), (i2,i3)
    tree = Phylo.read(io.StringIO("(diploid,(a,b)i1,(c,(d,e)i3)i2)i0;"), 'newick')
    names = sorted(c.name for c in tree.find_clades())
    original = Phylo.read(io.StringIO("(diploid,(a,b)i1,(c,(d,e)i3)i2)i0;"), 'newick')

    neighbours = list(nni_neighbours(tree))
    assert len(neighbours) == 6
    keys = {topology_key(n) for n in neighbours}
    assert len(keys) == 6
    assert topology_key(tree) not in keys
    for neighbour in neighbours:
        assert sorted(c.name for c in neighbour.find_clades()) == names
        assert neighbour.root.name == 'i0'
        assert 'diploid' in [c.name for c in neighbour.root.clades]
        assert all(len(c.clades) in (0, 2) or c is neighbour.root for c in neighbour.find_clades())
    assert topology_key(tree) == topology_key(original)  # input unchanged

    # the same moves on the output shape (unnamed root holding the normal and the MRCA)
    output_shape = Phylo.read(io.StringIO("(diploid,((a,b)i1,(c,(d,e)i3)i2)i0);"), 'newick')
    assert len(list(nni_neighbours(output_shape))) == 6


def test_moves_on_small_trees():
    assert list(nni_neighbours(read_tree("(diploid,(a,b)i1);"))) == []
    neighbours = list(nni_neighbours(read_tree("(diploid,((a,b)i2,c)i1);")))
    assert {topology_key(n) for n in neighbours} == {
        topology_key(read_tree("(diploid,((c,b)i2,a)i1);")),
        topology_key(read_tree("(diploid,((a,c)i2,b)i1);"))}


def test_topology_key_ignores_internal_names_and_child_order():
    a = Phylo.read(io.StringIO("(diploid,(a,b)i1,(c,(d,e)i3)i2)i0;"), 'newick')
    b = Phylo.read(io.StringIO("(((e,d)x,c)y,diploid,(b,a)z)w;"), 'newick')
    c = Phylo.read(io.StringIO("(diploid,(a,c)i1,(b,(d,e)i3)i2)i0;"), 'newick')
    assert topology_key(a) == topology_key(b)
    assert topology_key(a) != topology_key(c)


def test_nni_improves_bad_tree_to_local_optimum(meso5):
    _, evaluate = meso5
    start = read_tree(BAD_TREE)
    best, trace = nni_search(start, evaluate)
    assert trace[0] == dict(step=0, tree_length=98, evaluated=0)
    assert trace[-1]['tree_length'] == 67 < trace[0]['tree_length']
    assert [t['step'] for t in trace] == list(range(len(trace)))
    assert all(x['tree_length'] > y['tree_length'] for x, y in zip(trace, trace[1:]))
    assert evaluate(best) == trace[-1]['tree_length']
    assert topology_key(best) == topology_key(read_tree(open(DATA / 'meso5.new').read()))
    assert all(evaluate(n) >= trace[-1]['tree_length'] for n in nni_neighbours(best))


def test_nni_max_iter(meso5):
    _, evaluate = meso5
    calls = []

    def counting(tree):
        calls.append(1)
        return evaluate(tree)

    best, trace = nni_search(read_tree(BAD_TREE), counting, max_iter=8)
    assert len(calls) == 1 + 8
    assert trace[-1]['evaluated'] <= 8
    assert trace[-1]['tree_length'] < 98
    assert evaluate(best) == trace[-1]['tree_length']
    with pytest.raises(ValueError):
        nni_search(read_tree(BAD_TREE), evaluate, max_iter=0)


def test_nni_in_main_independent_of_cores(meso5):
    input_df, _ = meso5
    results = {}
    for n_cores in (1, 4):
        trace = []
        _, _, start_tree, final_tree, output_df, events_df = medicc.main(
            input_df, medicc.io.read_fst(), input_tree=import_tree(BAD_TREE), n_cores=n_cores,
            reconstruct_events=True, nni=True, nni_trace=trace)
        summary = medicc.summarize_patient(final_tree, None, input_df.index.unique('sample_id'))
        assert summary['tree_length'] == trace[-1]['tree_length'] == 67
        assert len(events_df) == 67
        assert topology_key(start_tree) != topology_key(final_tree)
        branch_lengths = {c.name: c.branch_length for c in final_tree.find_clades()}
        results[n_cores] = (topology_key(final_tree), branch_lengths, output_df, trace)
    assert results[1][:2] == results[4][:2]
    assert results[1][3] == results[4][3]
    pd.testing.assert_frame_equal(results[1][2], results[4][2])


def test_nni_scores_equal_main_tree_length(meso5):
    """The NNI score of a tree is the tree length of a normal run on that tree."""
    input_df, evaluate = meso5
    for newick in (BAD_TREE, open(DATA / 'meso5.new').read()):
        final_tree = medicc.main(input_df, medicc.io.read_fst(), input_tree=import_tree(newick))[3]
        summary = medicc.summarize_patient(final_tree, None, input_df.index.unique('sample_id'))
        assert evaluate(read_tree(newick)) == evaluate(read_tree(newick), n_cores=3) == \
            summary['tree_length']


def test_nni_with_custom_fst():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    fst = medicc.io.read_fst(user_fst=ROOT / 'medicc' / 'objects' / 'wgd_asymm.fst')
    trace = []
    medicc.main(input_df, fst, input_tree=import_tree(BAD_TREE), nni=True, nni_max_iter=2,
                nni_trace=trace)
    assert trace[-1]['evaluated'] <= 2
    assert trace[-1]['tree_length'] <= trace[0]['tree_length']


def test_main_rejects_nni_without_ancestors(meso5):
    input_df, _ = meso5
    with pytest.raises(medicc.core.MEDICCError, match='ancestral reconstruction'):
        medicc.main(input_df, medicc.io.read_fst(), ancestral_reconstruction=False, nni=True)
    with pytest.raises(medicc.core.MEDICCError, match='at least 1'):
        medicc.main(input_df, medicc.io.read_fst(), nni=True, nni_max_iter=0)


def test_cli_nni_writes_trace(tmp_path):
    tree_file = tmp_path / 'bad.new'
    tree_file.write_text(BAD_TREE)
    out = tmp_path / 'nni'
    run = run_medicc2(DATA / 'meso5.tsv', out, '--tree', tree_file, '--nni-mode', '--plot', 'none',
                      '--no-plot-tree', '--prefix', 'case', '-v')
    assert run.returncode == 0, run.stderr
    trace = pd.read_csv(out / 'case_nni_trace.tsv', sep='\t')
    assert list(trace.columns) == ['step', 'tree_length', 'evaluated']
    assert trace.iloc[0].tolist() == [0, 98, 0]
    assert trace['tree_length'].iloc[-1] == 67
    assert run.stderr.count('INFO medicc.nni:\tNNI step') == len(trace) - 1
    summary = pd.read_csv(out / 'case_summary.tsv', sep='\t', header=None, index_col=0)[1]
    assert float(summary['tree_length']) == 67

    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'plain', '--tree', tree_file, '--plot', 'none',
                      '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    assert not (tmp_path / 'plain' / 'case_nni_trace.tsv').exists()


@pytest.mark.parametrize('flags, message', [
    (['--nni-mode', '--topology-only'], '--topology-only'),
    (['--nni-mode', '--nni-max-iter', '0'], 'at least 1'),
    (['--nni-max-iter', '5'], 'requires --nni-mode')])
def test_cli_rejects_invalid_nni_options(tmp_path, flags, message):
    out = tmp_path / 'invalid'
    run = run_medicc2(DATA / 'meso5.tsv', out, *flags)
    assert run.returncode == 2
    assert message in run.stderr
    assert not out.exists()
