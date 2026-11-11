"""Tests of the tree inference methods FastME and neighbour joining (medicc.core.infer_tree_topology)."""
import io
import os
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

import medicc
from medicc.bootstrap import compare_trees
from medicc.conftest import find_test_fastme
from medicc.core import MEDICCError

ROOT = pathlib.Path(__file__).parent.parent.absolute()
DATA = pathlib.Path(__file__).parent / 'test_data'
MISSING = '/nonexistent/fastme'

# additive distances of the unrooted tree ((a:3,b:4):100000000,(c:7,(d:1,e:2):6):5,diploid:11)
ADDITIVE_LABELS = ['diploid', 'a', 'b', 'c', 'd', 'e']
ADDITIVE_TREE = "(diploid:11,((a:3,b:4)x:100000000,(c:7,(d:1,e:2)y:6)z:5)w:0);"


def _newick(text):
    return io.StringIO(text)


def additive_distances():
    """Distance matrix of ADDITIVE_TREE in the order of ADDITIVE_LABELS."""
    tree = Phylo.read(_newick(ADDITIVE_TREE), 'newick')
    return np.array([[tree.distance(x, y) for y in ADDITIVE_LABELS] for x in ADDITIVE_LABELS])


def run_medicc2(*args):
    return subprocess.run([sys.executable, 'medicc2', *map(str, args)], capture_output=True, text=True,
                          cwd=ROOT, timeout=600, env={**os.environ, 'MPLBACKEND': 'Agg'})


def leaf_sets(tree):
    """Set of leaf sets below each internal node of a rooted tree (topology without names)."""
    return {frozenset(c.name for c in clade.get_terminals()) for clade in tree.get_nonterminals()}


@pytest.mark.fastme
def test_fastme_main_tree_has_normal_as_outgroup():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    labels, _, _, tree, _, _ = medicc.main(input_df, medicc.io.read_fst(),
                                           fastme=find_test_fastme())
    assert sorted(c.name for c in tree.get_terminals()) == sorted(labels)
    assert 'diploid' in [c.name for c in tree.root.clades]
    assert all(len(c.clades) == 2 for c in tree.get_nonterminals())


@pytest.mark.fastme
def test_fastme_long_sample_names_round_trip():
    distances = additive_distances()
    long_labels = ['diploid'] + [name + '_' + 'x' * 70 for name in ADDITIVE_LABELS[1:]]
    assert all(len(name.encode('utf-8')) >= 64 for name in long_labels[1:])
    short = medicc.core.infer_tree_topology(distances, ADDITIVE_LABELS, 'diploid',
                                            fastme=find_test_fastme())
    long = medicc.core.infer_tree_topology(distances, long_labels, 'diploid',
                                           fastme=find_test_fastme())
    assert sorted(c.name for c in long.get_terminals()) == sorted(long_labels)
    for clade in long.get_terminals():
        clade.name = clade.name.split('_')[0]
    assert compare_trees(short, long)


@pytest.mark.fastme
@pytest.mark.parametrize('name', ['a b', 'a\tb', 'a(b', 'a)b', 'a,b', 'a:b', 'a;b', 'a[b', 'a]b', "a'b"])
def test_fastme_accepts_names_with_whitespace_or_newick_punctuation(name):
    labels = ADDITIVE_LABELS[:1] + [name] + ADDITIVE_LABELS[2:]
    tree = medicc.core.infer_tree_topology(additive_distances(), labels, 'diploid',
                                           fastme=find_test_fastme())
    assert sorted(c.name for c in tree.get_terminals()) == sorted(labels)
    reference = medicc.core.infer_tree_topology(additive_distances(), ADDITIVE_LABELS, 'diploid',
                                                fastme=find_test_fastme())
    for clade in tree.get_terminals():
        if clade.name == name:
            clade.name = 'a'
    assert compare_trees(reference, tree)


def test_unusable_fastme_falls_back_to_neighbour_joining(tmp_path, caplog):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    fst = medicc.io.read_fst()
    nj_tree = medicc.main(input_df, fst, tree_method='nj')[3]
    not_executable = tmp_path / 'not_executable'
    not_executable.write_text('#!/bin/sh\n')
    failing = tmp_path / 'failing'
    failing.write_text('#!/bin/sh\nexit 1\n')
    failing.chmod(0o755)
    for fastme, message in ((MISSING, 'FastME was not found'),
                            (not_executable, 'FastME was not found'),
                            (failing, 'FastME could not be used')):
        caplog.clear()
        with caplog.at_level('WARNING', logger='medicc'):
            tree = medicc.main(input_df, fst, fastme=fastme)[3]
        assert message in caplog.text
        assert 'using neighbour joining instead' in caplog.text
        assert compare_trees(tree, nj_tree)
        assert tree.total_branch_length() == nj_tree.total_branch_length()


def test_bootstrap_warns_once_without_fastme(caplog):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    with caplog.at_level('WARNING', logger='medicc'):
        medicc.bootstrap.run_bootstrap(input_df, N_bootstrap=5, show_progress=False, fastme=MISSING)
    assert caplog.text.count('FastME was not found') == 1


def test_fewer_than_four_samples_warn_and_use_neighbour_joining(caplog):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    three = input_df.loc[['diploid', 'leaf_3', 'leaf_4']]
    with caplog.at_level('WARNING', logger='medicc'):
        medicc.main(three, medicc.io.read_fst())
    assert 'FastME needs at least 4 samples (3 given)' in caplog.text


def test_cli_without_fastme_warns_and_uses_neighbour_joining(tmp_path):
    flags = ('--plot', 'none', '--no-plot-tree', '--prefix', 'case')
    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'missing', '--fastme', MISSING, *flags)
    assert run.returncode == 0, run.stderr
    log = run.stdout + run.stderr
    assert 'WARNING' in log and 'using neighbour joining instead' in log
    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'nj', '--tree-method', 'nj', *flags)
    assert run.returncode == 0, run.stderr
    for name in ('case_final_tree.new', 'case_branch_lengths.tsv', 'case_summary.tsv'):
        assert (tmp_path / 'missing' / name).read_text() == (tmp_path / 'nj' / name).read_text()


def test_tree_input_does_not_need_fastme(tmp_path):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    tree = medicc.io.import_tree(DATA / 'meso5.new')
    final_tree = medicc.main(input_df, medicc.io.read_fst(), input_tree=tree, fastme=MISSING)[3]
    assert final_tree.total_branch_length() == 67

    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'out', '--tree', DATA / 'meso5.new',
                      '--fastme', MISSING, '--plot', 'none', '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    summary = pd.read_csv(tmp_path / 'out/case_summary.tsv', sep='\t', header=None, index_col=0)[1]
    assert float(summary['tree_length']) == 67


def test_two_and_three_samples_do_not_call_fastme():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    fst = medicc.io.read_fst()
    two = input_df.loc[['diploid', 'leaf_0']]
    _, distances, _, tree, _, _ = medicc.main(two, fst, fastme=MISSING)
    assert tree.total_branch_length() == distances.loc['diploid', 'leaf_0'] == 19

    three = input_df.loc[['diploid', 'leaf_3', 'leaf_4']]
    _, distances, _, tree, _, _ = medicc.main(three, fst, fastme=MISSING)
    assert sorted(c.name for c in tree.get_terminals()) == ['diploid', 'leaf_3', 'leaf_4']
    assert 'diploid' in [c.name for c in tree.root.clades]
    fixed = medicc.io.import_tree(_newick("(diploid,(leaf_3,leaf_4)internal_1);"))
    assert tree.total_branch_length() == medicc.main(three, fst, input_tree=fixed)[3].total_branch_length()
    # the only unrooted topology, built by neighbour joining (three-point branch lengths)
    d = distances.values
    star = medicc.core.infer_tree_topology(d, list(distances.index), 'diploid', fastme=MISSING)
    lengths = {c.name: c.branch_length for c in star.root.clades}
    assert lengths['diploid'] == (d[0, 1] + d[0, 2] - d[1, 2]) / 2
    assert lengths['leaf_3'] == (d[0, 1] + d[1, 2] - d[0, 2]) / 2


@pytest.mark.fastme
def test_cli_two_and_three_samples(tmp_path):
    input_df = pd.read_csv(DATA / 'meso5.tsv', sep='\t')
    for samples in (['diploid', 'leaf_0'], ['diploid', 'leaf_3', 'leaf_4']):
        name = f"n{len(samples)}"
        input_df[input_df['sample_id'].isin(samples)].to_csv(tmp_path / f'{name}.tsv', sep='\t',
                                                             index=False)
        run = run_medicc2(tmp_path / f'{name}.tsv', tmp_path / name, '--events', '--plot', 'none',
                          '--no-plot-tree')
        assert run.returncode == 0, run.stderr
        tree = Phylo.read(tmp_path / name / f'{name}_final_tree.new', 'newick')
        assert sorted(c.name for c in tree.get_terminals()) == samples


@pytest.mark.fastme
def test_bootstrap_uses_the_given_fastme():
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    trees_df, _ = medicc.bootstrap.run_bootstrap(input_df, N_bootstrap=5, seed=1,
                                                 show_progress=False, fastme=find_test_fastme())
    assert trees_df['count'].sum() == 5
    for tree in trees_df['tree']:
        assert sorted(c.name for c in tree.get_terminals()) == \
            sorted(input_df.index.get_level_values('sample_id').unique())


def test_invalid_tree_method_is_rejected(tmp_path):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    with pytest.raises(MEDICCError, match="Unknown tree method 'upgma'"):
        medicc.main(input_df, medicc.io.read_fst(), tree_method='upgma')
    with pytest.raises(MEDICCError, match='Unknown tree method'):
        medicc.core.infer_tree_topology(additive_distances(), ADDITIVE_LABELS, 'diploid',
                                        method='upgma')
    with pytest.raises(MEDICCError, match='Unknown tree method'):
        medicc.bootstrap.run_bootstrap(input_df, N_bootstrap=5, show_progress=False,
                                       tree_method='upgma')
    out = tmp_path / 'out'
    run = run_medicc2(DATA / 'meso5.tsv', out, '--tree-method', 'upgma')
    assert run.returncode == 2
    assert "invalid choice: 'upgma'" in run.stderr
    assert not out.exists()


def test_nj_method_does_not_need_fastme(tmp_path):
    input_df = medicc.io.read_and_parse_input_data(DATA / 'meso5.tsv')
    labels, distances, _, tree, _, _ = medicc.main(input_df, medicc.io.read_fst(),
                                                   tree_method='nj', fastme=MISSING)
    expected = medicc.nj.NeighbourJoining(distances.values, distances.index).tree
    expected.root_with_outgroup('diploid')
    tree.root_with_outgroup('diploid')
    assert leaf_sets(tree) == leaf_sets(expected)
    trees_df, _ = medicc.bootstrap.run_bootstrap(input_df, N_bootstrap=5, seed=1, show_progress=False,
                                                 tree_method='nj', fastme=MISSING)
    assert trees_df['count'].sum() == 5

    run = run_medicc2(DATA / 'meso5.tsv', tmp_path / 'out', '--tree-method', 'nj', '--fastme', MISSING,
                      '--bootstrap-nr', '5', '--plot', 'none', '--no-plot-tree', '--prefix', 'case')
    assert run.returncode == 0, run.stderr
    assert (tmp_path / 'out/case_support_tree.new').is_file()


@pytest.mark.fastme
@pytest.mark.parametrize('method', ['fastme', 'nj'])
def test_both_methods_recover_an_additive_tree_rooted_on_the_normal(method):
    tree = medicc.core.infer_tree_topology(additive_distances(), ADDITIVE_LABELS, 'diploid',
                                           method=method, fastme=find_test_fastme())
    assert sorted(c.name for c in tree.get_terminals()) == sorted(ADDITIVE_LABELS)
    assert sorted(c.name for c in tree.get_nonterminals()) == [f'internal{i}' for i in range(1, 5)]
    assert len(tree.root.clades) == 3 and 'diploid' in [c.name for c in tree.root.clades]
    tree.root_with_outgroup('diploid')
    expected = Phylo.read(_newick(ADDITIVE_TREE), 'newick')
    expected.root_with_outgroup('diploid')
    assert leaf_sets(tree) == leaf_sets(expected)
    assert tree.total_branch_length() == pytest.approx(expected.total_branch_length())


# OV03-04 final trees (normal OV03-04_diploid): neighbour joining as in MEDICC2 before FastME
# (tree length 727) and FastME (tree length 724, other topology)
OV03_04_NJ_TREE = (
    '((((((((OV03-04_Om_S03:23,(OV03-04_Om_S02:17,(OV03-04_Om_S06:27,'
    'OV03-04_Om_S09:17)internal14:8)internal16:6)internal17:12,OV03-04_Om_S07:22)internal15:18,'
    'OV03-04_Om_S15:16)internal13:16,OV03-04_Om_S05:35)internal12:14,(OV03-04_Bl_B03:37,'
    '(OV03-04_VV_B01:22,(OV03-04_Bl_B01:25,'
    'OV03-04_VV_B02:18)internal7:12)internal8:20)internal10:29)internal11:22,((OV03-04_Om_S08:16,'
    'OV03-04_Om_S12:26)internal3:22,(OV03-04_Om_S01:28,(OV03-04_Om_S11:30,'
    'OV03-04_Om_S13:25)internal4:12)internal5:16)internal6:8)internal9:37,(OV03-04_Om_S04:16,'
    'OV03-04_Om_S14:12)internal1:15)internal2:48,OV03-04_diploid:0):0;')
OV03_04_FASTME_TREE = (
    '(((((((OV03-04_VV_B02:18,OV03-04_Bl_B01:25)internal17:12,OV03-04_VV_B01:22)internal16:19,'
    'OV03-04_Bl_B03:38)internal15:24,(((OV03-04_Om_S02:13,(OV03-04_Om_S03:24,'
    'OV03-04_Om_S15:18)internal1:10)internal2:13,OV03-04_Om_S07:23)internal3:6,'
    '(OV03-04_Om_S06:25,OV03-04_Om_S09:19)internal4:12)internal5:22)internal14:26,'
    '((OV03-04_Om_S08:14,OV03-04_Om_S12:28)internal6:24,((OV03-04_Om_S01:30,'
    'OV03-04_Om_S11:29)internal7:16,(OV03-04_Om_S05:38,'
    'OV03-04_Om_S13:18)internal8:11)internal9:12)internal10:8)internal13:37,(OV03-04_Om_S04:16,'
    'OV03-04_Om_S14:12)internal11:15)internal12:47,OV03-04_diploid:0):0;')


def clade_values(tree):
    """Rooted tree as {leaf set: (node name, branch length)}."""
    return {frozenset(c.name for c in clade.get_terminals()): (clade.name, clade.branch_length)
            for clade in tree.find_clades()}


@pytest.fixture(scope='module')
def ov03_04():
    return medicc.io.read_and_parse_input_data(
        ROOT / 'examples/OV03-04/OV03-04_descr.txt', normal_name='OV03-04_diploid', input_type='fasta')


def test_nj_method_reproduces_previous_trees(ov03_04):
    final_tree = medicc.main(ov03_04, medicc.io.read_fst(), normal_name='OV03-04_diploid',
                             n_cores=2, tree_method='nj', fastme=MISSING)[3]
    assert final_tree.total_branch_length() == 727
    assert clade_values(final_tree) == clade_values(Phylo.read(_newick(OV03_04_NJ_TREE), 'newick'))


@pytest.mark.fastme
def test_fastme_method_tree_of_ov03_04(ov03_04):
    final_tree = medicc.main(ov03_04, medicc.io.read_fst(), normal_name='OV03-04_diploid',
                             n_cores=2, fastme=find_test_fastme())[3]
    assert final_tree.total_branch_length() == 724
    assert clade_values(final_tree) == clade_values(Phylo.read(_newick(OV03_04_FASTME_TREE), 'newick'))
