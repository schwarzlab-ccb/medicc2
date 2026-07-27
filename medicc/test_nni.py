import fstlib
import numpy as np
import pytest
from io import StringIO

import medicc
import Bio.Phylo as Phylo

def _fresh_tree(tree_newick):
  return Phylo.read(StringIO(tree_newick), "newick")


def test_apply_nni_swap():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")
    u_name = "n3"
    v_name = "n1"
    A_name = "n2"
    swap_target_name = "B"

    target_tree_newick = "((((A,(C,D)n2)n1,B)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;"

    nni_tree_after_swap = medicc.nni._apply_nni_swap(source_tree, u_name, v_name, A_name, swap_target_name)
    target_tree = _fresh_tree(target_tree_newick)

    nni_tree_newick_after_swap = medicc.tree_hash.get_canonical_newick(nni_tree_after_swap)
    target_tree_newick_after_swap = medicc.tree_hash.get_canonical_newick(target_tree)

    assert nni_tree_newick_after_swap == target_tree_newick_after_swap

def test_enumerate_moves():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")
    moves = medicc.nni._enumerate_moves(source_tree)

    expected_moves = [
        medicc.nni._NNIMove(u_name="mrca", v_name="n3", A_name="n6", swap_target_name="n1"),
        medicc.nni._NNIMove(u_name="mrca", v_name="n3", A_name="n6", swap_target_name="n2"),
        medicc.nni._NNIMove(u_name="mrca", v_name="n6", A_name="n3", swap_target_name="n4"),
        medicc.nni._NNIMove(u_name="mrca", v_name="n6", A_name="n3", swap_target_name="n5"),
        medicc.nni._NNIMove(u_name="n3", v_name="n1", A_name="n2", swap_target_name="A"),
        medicc.nni._NNIMove(u_name="n3", v_name="n1", A_name="n2", swap_target_name="B"),
        medicc.nni._NNIMove(u_name="n3", v_name="n2", A_name="n1", swap_target_name="C"),
        medicc.nni._NNIMove(u_name="n3", v_name="n2", A_name="n1", swap_target_name="D"),
        medicc.nni._NNIMove(u_name="n6", v_name="n4", A_name="n5", swap_target_name="E"),
        medicc.nni._NNIMove(u_name="n6", v_name="n4", A_name="n5", swap_target_name="F"),
        medicc.nni._NNIMove(u_name="n6", v_name="n5", A_name="n4", swap_target_name="G"),
        medicc.nni._NNIMove(u_name="n6", v_name="n5", A_name="n4", swap_target_name="H"),
    ]

    assert len(moves) == len(expected_moves)
    for move in expected_moves:
        assert any(
            m.u_name == move.u_name and
            m.v_name == move.v_name and
            m.A_name == move.A_name and
            m.swap_target_name == move.swap_target_name
            for m in moves
        )

def test_nni_neighbors_count_matches_enumerate_moves():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")

    moves = medicc.nni._enumerate_moves(source_tree)
    neighbors = list(medicc.nni.nni_neighbors(source_tree))

    assert len(neighbors) == len(moves)


def test_nni_neighbors_matches_apply_nni_swap():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")

    for new_tree, move in medicc.nni.nni_neighbors(source_tree):
        expected_tree = medicc.nni._apply_nni_swap(
            source_tree, move.u_name, move.v_name, move.A_name, move.swap_target_name)

        new_tree_newick = medicc.tree_hash.get_canonical_newick(new_tree)
        expected_tree_newick = medicc.tree_hash.get_canonical_newick(expected_tree)

        assert new_tree_newick == expected_tree_newick


def test_nni_neighbors_does_not_mutate_input_tree():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")
    original_newick = medicc.tree_hash.get_canonical_newick(source_tree)

    list(medicc.nni.nni_neighbors(source_tree))

    assert medicc.tree_hash.get_canonical_newick(source_tree) == original_newick


def test_nni_neighbors_produces_distinct_topologies():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")

    hashes = [medicc.tree_hash.get_topology_hash(new_tree)
             for new_tree, _ in medicc.nni.nni_neighbors(source_tree)]

    assert len(hashes) == len(set(hashes))


def test_nni_neighbors_specific_swap():
    source_tree = _fresh_tree("((((A,B)n1,(C,D)n2)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")
    target_tree = _fresh_tree("((((A,(C,D)n2)n1,B)n3,((E,F)n4,(G,H)n5)n6)mrca)diploid;")
    target_newick = medicc.tree_hash.get_canonical_newick(target_tree)

    matches = [
        new_tree for new_tree, move in medicc.nni.nni_neighbors(source_tree)
        if move.u_name == "n3" and move.v_name == "n1"
        and move.A_name == "n2" and move.swap_target_name == "B"
    ]

    assert len(matches) == 1
    assert medicc.tree_hash.get_canonical_newick(matches[0]) == target_newick
