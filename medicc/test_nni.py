import copy
import logging
import pathlib

import fstlib
import numpy as np
import pytest
from io import StringIO

import medicc
import Bio.Phylo as Phylo

def _fresh_tree(tree_newick):
  return Phylo.read(StringIO(tree_newick), "newick")

FST = medicc.io.read_fst()

def _fsa(seq):
    return fstlib.factory.from_string(seq, isymbols=FST.input_symbols(), osymbols=FST.output_symbols(), arc_type=FST.arc_type())


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


# --- _reshape_nj_for_search / _wrap_tree_for_output ---

def test_reshape_nj_for_search_produces_search_shape():
    nj_tree = _fresh_tree("(diploid,(A,B)mrca)root;")

    reshaped = medicc.core._reshape_nj_for_search(nj_tree, "diploid")

    assert reshaped.root.name == "diploid"
    assert len(reshaped.root.clades) == 1
    assert reshaped.root.clades[0].name == "mrca"
    assert [l.name for l in reshaped.get_terminals()] == ["A", "B"]


def test_reshape_nj_for_search_does_not_mutate_input():
    nj_tree = _fresh_tree("(diploid,(A,B)mrca)root;")
    original_newick = medicc.tree_hash.get_canonical_newick(nj_tree)

    medicc.core._reshape_nj_for_search(nj_tree, "diploid")

    assert medicc.tree_hash.get_canonical_newick(nj_tree) == original_newick


def test_wrap_tree_for_output_restores_standard_shape():
    normal_name = "diploid"
    samples = {
        "diploid": _fsa("11111X1111"),
        "A": _fsa("22111X1111"),
        "B": _fsa("21211X1111"),
    }
    nj_tree = _fresh_tree("(diploid,(A,B)mrca)root;")
    search_tree = medicc.core._reshape_nj_for_search(nj_tree, normal_name)
    ancestors, _ = medicc.ancestors.reconstruct_ancestors(
        search_tree, samples, FST, FST, normal_name, upper_cache=True)

    wrapped = medicc.core._wrap_tree_for_output(search_tree, FST, ancestors, normal_name)

    assert wrapped.root.name is None
    assert len(wrapped.root.clades) == 2
    diploid_clade = next(c for c in wrapped.root.clades if c.name == normal_name)
    mrca_clade = next(c for c in wrapped.root.clades if c.name != normal_name)
    assert diploid_clade.clades == []
    assert [l.name for l in mrca_clade.get_terminals()] == ["A", "B"]


def test_wrap_tree_for_output_recomputes_branch_lengths_before_reroot():
    # _wrap_tree_for_output must score the mrca edge under `fst` while it is
    # still attached to the search-shape (diploid-named) root; after the
    # re-root its parent becomes the new unnamed root, which
    # update_branch_lengths would then skip over instead of scoring.
    normal_name = "diploid"
    samples = {
        "diploid": _fsa("11111X1111"),
        "A": _fsa("22111X1111"),
        "B": _fsa("21211X1111"),
    }
    nj_tree = _fresh_tree("(diploid,(A,B)mrca)root;")
    search_tree = medicc.core._reshape_nj_for_search(nj_tree, normal_name)
    ancestors, _ = medicc.ancestors.reconstruct_ancestors(
        search_tree, samples, FST, FST, normal_name, upper_cache=True)

    reference_tree = copy.deepcopy(search_tree)
    medicc.core.update_branch_lengths(reference_tree, FST, ancestors, normal_name)
    expected_mrca_branch_length = reference_tree.root.clades[0].branch_length
    assert expected_mrca_branch_length > 0  # otherwise this test can't discriminate

    wrapped = medicc.core._wrap_tree_for_output(search_tree, FST, ancestors, normal_name)
    mrca_clade = next(c for c in wrapped.root.clades if c.name != normal_name)

    assert mrca_clade.branch_length == expected_mrca_branch_length


# --- nni_mode hill-climb state machine ---

NORMAL = "diploid"
TREE4 = "(((A,B)n1,(C,D)n2)mrca)diploid;"


def _run_nni_mode(profiles, max_iter=5000):
    tree = _fresh_tree(TREE4)
    samples = {"diploid": _fsa("1111X1111")}
    samples.update({name: _fsa(profile) for name, profile in profiles.items()})
    reshaped = medicc.core._reshape_nj_for_search(tree, NORMAL)
    return medicc.core.nni_mode(reshaped, samples, FST, FST, normal_name=NORMAL, nni_max_iter=max_iter)


def test_nni_mode_terminates_immediately_when_no_moves_exist():
    # a cherry has no internal-internal edge, so _enumerate_moves is always empty
    tree = _fresh_tree("(A,B)diploid;")
    reshaped = medicc.core._reshape_nj_for_search(tree, NORMAL)
    samples = {"diploid": _fsa("11X11"), "A": _fsa("21X11"), "B": _fsa("11X21")}

    result = medicc.core.nni_mode(reshaped, samples, FST, FST, normal_name=NORMAL, nni_max_iter=5000)

    assert len(result["trace"]) == 1  # only the seeded initial score, no sweep ran
    assert len(result["best_trees"]) == 1
    assert medicc.tree_hash.get_topology_hash(result["best_trees"][0]) == \
        medicc.tree_hash.get_topology_hash(reshaped)


def test_nni_mode_converges_to_unique_optimum_via_strict_improvements():
    crossed = {"A": "2111X1111", "C": "2111X1111", "B": "1111X2111", "D": "1111X2111"}

    result = _run_nni_mode(crossed)

    assert result["trace"] == [4.0, 3.0, 2.0]
    assert result["best_score"] == 2.0
    assert len(result["best_trees"]) == 1


def test_nni_mode_trace_is_strictly_decreasing():
    crossed = {"A": "2111X1111", "C": "2111X1111", "B": "1111X2111", "D": "1111X2111"}

    trace = _run_nni_mode(crossed)["trace"]

    assert all(trace[i] > trace[i + 1] for i in range(len(trace) - 1))


def test_nni_mode_plateau_expands_across_sweeps_to_distinct_topologies():
    identical = {"A": "1111X1111", "B": "1111X1111", "C": "1111X1111", "D": "1111X1111"}

    result = _run_nni_mode(identical, max_iter=5000)

    assert result["best_score"] == 0.0
    assert len(result["best_trees"]) == 15
    hashes = [medicc.tree_hash.get_topology_hash(t) for t in result["best_trees"]]
    assert len(hashes) == len(set(hashes))


def test_nni_mode_hard_cap_warns_and_returns_partial_plateau(caplog):
    identical = {"A": "1111X1111", "B": "1111X1111", "C": "1111X1111", "D": "1111X1111"}

    with caplog.at_level(logging.WARNING, logger="medicc.core"):
        result = _run_nni_mode(identical, max_iter=9)

    assert len(result["best_trees"]) == 9  # less than the natural 15-tree plateau
    assert any("hard cap hit" in record.message for record in caplog.records)


def test_nni_mode_hard_cap_can_interrupt_mid_sweep():
    # Regression test for the deliberate 2026-07-28 change (commit c122827):
    # nni_max_iter is checked once per frontier tree, not once per whole sweep,
    # so it can stop partway through a sweep's frontier rather than only
    # between sweeps. Several small caps landing in the same "gap" must give
    # the identical truncated result; crossing the next threshold must not.
    identical = {"A": "1111X1111", "B": "1111X1111", "C": "1111X1111", "D": "1111X1111"}

    n_at_9 = len(_run_nni_mode(identical, max_iter=9)["best_trees"])
    n_at_12 = len(_run_nni_mode(identical, max_iter=12)["best_trees"])
    n_at_13 = len(_run_nni_mode(identical, max_iter=13)["best_trees"])

    assert n_at_9 == n_at_12 == 9
    assert n_at_13 == 11


# --- main()'s nni_export_all_topology wiring ---

EXAMPLE_TSV = str(pathlib.Path(__file__).resolve().parent.parent / "examples" / "simple_example" / "simple_example.tsv")


def _example_input_df():
    return medicc.io.read_and_parse_input_data(
        filename=EXAMPLE_TSV,
        normal_name=NORMAL,
        input_type="TSV",
        separator="X",
        allele_columns=["cn_a", "cn_b"],
        total_copy_numbers=False,
        maxcn=8)


def test_main_nni_single_topology_matches_first_export_all_topology():
    input_df = _example_input_df()

    single_result = medicc.core.main(
        input_df=input_df, asymm_upper_fst=FST, asymm_lower_fst=FST, normal_name=NORMAL,
        chr_separator="X", nni_mode_flag=True, nni_export_all_topology=False,
        reconstruct_events=True, nni_max_iter=2000)
    all_result = medicc.core.main(
        input_df=input_df, asymm_upper_fst=FST, asymm_lower_fst=FST, normal_name=NORMAL,
        chr_separator="X", nni_mode_flag=True, nni_export_all_topology=True,
        reconstruct_events=True, nni_max_iter=2000)

    single_tree, single_output_df = single_result[3], single_result[4]
    tree_l, output_df_l = all_result[3], all_result[4]

    assert single_tree.format("newick") == tree_l[0].format("newick")
    assert single_output_df.equals(output_df_l[0])


def test_main_nni_export_all_topologies_are_distinct():
    input_df = _example_input_df()

    result = medicc.core.main(
        input_df=input_df, asymm_upper_fst=FST, asymm_lower_fst=FST, normal_name=NORMAL,
        chr_separator="X", nni_mode_flag=True, nni_export_all_topology=True,
        reconstruct_events=True, nni_max_iter=2000)
    tree_l, output_df_l, events_df_l = result[3], result[4], result[5]

    assert len(tree_l) > 1  # confirms this fixture genuinely plateaus
    newicks = [t.format("newick") for t in tree_l]
    assert len(set(newicks)) == len(newicks)
    assert len({id(df) for df in output_df_l}) == len(output_df_l)
    assert len({id(df) for df in events_df_l}) == len(events_df_l)


def test_main_nni_export_all_without_reconstruct_events_does_not_crash():
    # Regression test: reconstruct_events defaults to False, which previously
    # left events_df_l undefined and crashed the final return with a NameError.
    input_df = _example_input_df()

    result = medicc.core.main(
        input_df=input_df, asymm_upper_fst=FST, asymm_lower_fst=FST, normal_name=NORMAL,
        chr_separator="X", nni_mode_flag=True, nni_export_all_topology=True,
        nni_max_iter=2000)

    assert result[5] is None


def test_main_plain_path_return_arity_is_backward_compatible():
    # Regression test: main()'s return grew from 6 to 8 values for NNI mode;
    # the plain (non-NNI) path must keep returning exactly 6 so pre-existing
    # callers (the medicc2 CLI, bootstrap.py) keep unpacking correctly.
    input_df = _example_input_df()

    result = medicc.core.main(
        input_df=input_df, asymm_upper_fst=FST, asymm_lower_fst=FST, normal_name=NORMAL,
        chr_separator="X", nni_mode_flag=False)

    assert len(result) == 6
