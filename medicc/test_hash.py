import fstlib
import numpy as np
import pytest
from io import StringIO

import medicc
import Bio.Phylo as Phylo

def _fresh_tree(tree_newick):
  return Phylo.read(StringIO(tree_newick), "newick")

def test_identical_trees():
    newick_1 = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"
    identical_newick_1 = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"

    tree_1 = _fresh_tree(newick_1)
    tree_2 = _fresh_tree(identical_newick_1)

    tree_1_hash = medicc.tree_hash.get_topology_hash(tree_1)
    tree_2_hash = medicc.tree_hash.get_topology_hash(tree_2)

    assert tree_1_hash == tree_2_hash

def test_different_ancestor_name():
    newick_1 = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"
    newick_1_rename_branch_length = "(((A,B)I37,(C,D)I47)I87,((E,F)I57,(G,H)I67)I27)diploid;"

    tree_1 = _fresh_tree(newick_1)
    tree_2 = _fresh_tree(newick_1_rename_branch_length)

    tree_1_hash = medicc.tree_hash.get_topology_hash(tree_1)
    tree_2_hash = medicc.tree_hash.get_topology_hash(tree_2)

    assert tree_1_hash == tree_2_hash

def test_with_branch_length_difference():
    newick_1 = "(((A:1,B:2)I3:3,(C:4,D:5)I4:6)I1:7,((E:8,F:9)I5:10,(G:11,H:12)I6:13)I2:14)diploid;"
    newick_1_rename_branch_length = "(((A:10,B:20)I37:30,(C:40,D:50)I47:60)I87:70,((E:80,F:90)I57:100,(G:110,H:120)I67:130)I27:140)diploid;"

    tree_1 = _fresh_tree(newick_1)
    tree_2 = _fresh_tree(newick_1_rename_branch_length)

    tree_1_hash = medicc.tree_hash.get_topology_hash(tree_1)
    tree_2_hash = medicc.tree_hash.get_topology_hash(tree_2)

    assert tree_1_hash == tree_2_hash

def test_swapping_child_order():
    newick_1 = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"
    newick_1_swap_child_order = "(((B,A)I3,(D,C)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"

    tree_1 = _fresh_tree(newick_1)
    tree_2 = _fresh_tree(newick_1_swap_child_order)

    tree_1_hash = medicc.tree_hash.get_topology_hash(tree_1)
    tree_2_hash = medicc.tree_hash.get_topology_hash(tree_2)

    assert tree_1_hash == tree_2_hash

def test_different_trees():
    newick_1 = "(((A,B)I3,(C,D)I4)I1,((E,F)I5,(G,H)I6)I2)diploid;"
    newick_2 = "(A,(B,(C,(D,E)I3)I2)I1)Root;"

    tree_1 = _fresh_tree(newick_1)
    tree_2 = _fresh_tree(newick_2)

    tree_1_hash = medicc.tree_hash.get_topology_hash(tree_1)
    tree_2_hash = medicc.tree_hash.get_topology_hash(tree_2)

    assert tree_1_hash != tree_2_hash

def test_subtree_swap_changes_topology():
  # before: u's children are {A, v}, v's children are {B, C}
  before = "((A,(B,C)v)u,X)diploid;"
  # after an NNI-style swap of A and B: u's children become {B, v}, v's children become {A, C}
  after  = "((B,(A,C)v)u,X)diploid;"

  tree_1 = _fresh_tree(before)
  tree_2 = _fresh_tree(after)

  tree_1_hash = medicc.tree_hash.get_topology_hash(tree_1)
  tree_2_hash = medicc.tree_hash.get_topology_hash(tree_2)

  assert tree_1_hash != tree_2_hash


def test_canonical_newick_is_not_branch_length_invariant():
  a = _fresh_tree("(((A:1,B:2)I3:3,(C:4,D:5)I4:6)I1:7,X:8)diploid;")
  b = _fresh_tree("(((A:9,B:9)I3:9,(C:9,D:9)I4:9)I1:9,X:9)diploid;")
  assert medicc.tree_hash.get_canonical_newick(a) != medicc.tree_hash.get_canonical_newick(b)
  # but topology hash still agrees:
  assert medicc.tree_hash.get_topology_hash(a) == medicc.tree_hash.get_topology_hash(b)


