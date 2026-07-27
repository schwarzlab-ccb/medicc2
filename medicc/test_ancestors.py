from dataclasses import dataclass

from setuptools.command.bdist_egg import iter_symbols

import fstlib
import numpy as np
import pytest
from io import StringIO

import medicc
import Bio.Phylo as Phylo

FST = medicc.io.read_fst()

def __extract_level_name(levels):
    level_name = {}
    for level in levels:
        level_name[level] = [i.name for i in levels[level]]
    return level_name

def _fsa(seq):
    return fstlib.factory.from_string(seq, isymbols=FST.input_symbols(), osymbols=FST.output_symbols(), arc_type=FST.arc_type())

def _decode_all(fsa_dict):
  return {name: medicc.tools.fsa_to_string(fsa) for name, fsa in fsa_dict.items()}

NORMAL = "diploid"
TREE_NEWICK = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"

def _fresh_tree():
  return Phylo.read(StringIO(TREE_NEWICK), "newick")

SAMPLES_DICT = {
  "diploid": _fsa("11111X1111"),
  "A": _fsa("11111X1111"), "B": _fsa("21111X1111"),
  "C": _fsa("11111X1111"), "D": _fsa("11101X1111"),
  "E": _fsa("22111X1111"), "F": _fsa("11111X2111"),
  "G": _fsa("11111X1111"), "H": _fsa("11121X1111"),
}

def test_group_nodes_by_depth_caterpillar_tree():
    normal_name = "Root"
    caterpillar_tree_newick = "(A,(B,(C,(D,E)I3)I2)I1)Root;"
    caterpillar_tree = Phylo.read(StringIO(caterpillar_tree_newick), "newick")

    clade_list = [clade for clade in caterpillar_tree.find_clades(order="preorder") if clade.name != normal_name]
    levels = medicc.ancestors._group_nodes_by_depth(clade_list, normal_name)
    expected_levels = {
        0: ['I1'],
        1: ['I2'],
        2: ['I3'],
    }
    levels_name = __extract_level_name(levels)

    assert levels_name == expected_levels

def test_group_nodes_by_depth_perfectly_balanced_tree():
    normal_name = "Root"
    perfectly_balanced_tree_newick = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)Root;"
    perfectly_balanced_tree = Phylo.read(StringIO(perfectly_balanced_tree_newick), "newick")

    clade_list = [clade for clade in perfectly_balanced_tree.find_clades(order="preorder") if clade.name != normal_name]
    levels = medicc.ancestors._group_nodes_by_depth(clade_list, normal_name)

    expected_levels = {
        0: ['I1', 'I2'],
        1: ['I3', 'I4', 'I5', 'I6'],
    }

    levels_name = __extract_level_name(levels)
    assert levels_name == expected_levels

def test_group_nodes_by_depth_irregular_binary_tree():
    normal_name = "Root"
    random_tree_newick = "((A,(B,C)I3)I1,((D,E)I4,(F,(G,H)I6)I5)I2)Root;"
    random_tree = Phylo.read(StringIO(random_tree_newick), "newick")

    clade_list = [clade for clade in random_tree.find_clades(order="preorder") if clade.name != normal_name]
    levels = medicc.ancestors._group_nodes_by_depth(clade_list, normal_name)

    expected_levels = {
        0: ['I1', 'I2'],
        1: ['I3', 'I4', 'I5'],
        2: ['I6'],
    }

    levels_name = __extract_level_name(levels)
    assert levels_name == expected_levels


@pytest.mark.parametrize("n_cores", [None, 1, 2, 4])
def test_reconstruct_ancestors_parallel_matches_serial(n_cores):
    baseline = medicc.reconstruct_ancestors(
        _fresh_tree(), SAMPLES_DICT, FST, FST, NORMAL, n_cores=None)
    result = medicc.reconstruct_ancestors(
        _fresh_tree(), SAMPLES_DICT, FST, FST, NORMAL, n_cores=n_cores)
    assert _decode_all(result) == _decode_all(baseline)

