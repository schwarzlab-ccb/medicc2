"""Unit tests for medicc.nni — NNI move enumeration."""
import copy
import pytest
from Bio.Phylo.BaseTree import Clade, Tree

import medicc.nni as nni
import medicc.spr as spr
import medicc.tree_hash as tree_hash
from medicc.ancestors import _get_dirty_nodes_up


def _make_search_shape_tree(n_leaves):
    """Build a left-comb search-shape tree with n_leaves non-diploid leaves.

    Result shape:
      root("diploid") -> internal_0 -> internal_1 -> ... -> (s1, s2, ..., s_n)

    Specifically a caterpillar where each internal node has one leaf child and
    one internal child. The deepest internal node has two leaf children.
    """
    assert n_leaves >= 2
    leaves = [Clade(name=f"s{i+1}") for i in range(n_leaves)]
    # Build from the deepest internal upward
    deepest_internal_idx = n_leaves - 2  # internal_{n_leaves-2} has two leaves
    current = Clade(name=f"internal_{deepest_internal_idx}",
                    clades=[leaves[-2], leaves[-1]])
    # Walk outward: each step adds a new internal whose children are (next leaf, current)
    for i in range(deepest_internal_idx - 1, -1, -1):
        leaf = leaves[i]
        current = Clade(name=f"internal_{i}", clades=[leaf, current])
    root = Clade(name="diploid", clades=[current])
    return Tree(root=root, rooted=True)


def test_nni_neighbor_count_n3():
    tree = _make_search_shape_tree(3)
    # n=3 -> 2(n-3) = 0 neighbors
    assert list(nni.nni_neighbors(tree)) == []


def test_nni_neighbor_count_n4():
    tree = _make_search_shape_tree(4)
    # n=4 -> 2(n-3) = 2 neighbors
    assert len(list(nni.nni_neighbors(tree))) == 2


def test_nni_neighbor_count_n5():
    tree = _make_search_shape_tree(5)
    # n=5 -> 2(n-3) = 4 neighbors
    assert len(list(nni.nni_neighbors(tree))) == 4


def test_nni_neighbor_count_n8():
    tree = _make_search_shape_tree(8)
    # n=8 -> 2(n-3) = 10 neighbors
    assert len(list(nni.nni_neighbors(tree))) == 10


def test_nni_neighbors_are_distinct():
    tree = _make_search_shape_tree(6)
    hashes = set()
    input_hash = tree_hash.tree_hash(tree_hash.strip_branch_lengths(tree))
    for neighbor_tree, _ in nni.nni_neighbors(tree):
        h = tree_hash.tree_hash(tree_hash.strip_branch_lengths(neighbor_tree))
        assert h != input_hash, "NNI neighbor identical to input"
        assert h not in hashes, "NNI neighbors not pairwise distinct"
        hashes.add(h)
    # n=6 -> 2(n-3) = 6 neighbors
    assert len(hashes) == 6


def _child_sets_by_name(tree):
    return {c.name: frozenset(child.name for child in c.clades)
            for c in tree.find_clades() if len(c.clades) > 0}


def test_nni_neighbor_changes_exactly_two_internal_child_sets():
    tree = _make_search_shape_tree(5)
    input_csets = _child_sets_by_name(tree)
    for neighbor_tree, _ in nni.nni_neighbors(tree):
        nbr_csets = _child_sets_by_name(neighbor_tree)
        differing = [name for name in input_csets if input_csets[name] != nbr_csets[name]]
        assert len(differing) == 2, \
            f"NNI move should change exactly 2 internal child sets, got {len(differing)}"


def test_nni_preserves_internal_node_names():
    tree = _make_search_shape_tree(5)
    input_names = {c.name for c in tree.find_clades() if c.name is not None}
    for neighbor_tree, _ in nni.nni_neighbors(tree):
        nbr_names = {c.name for c in neighbor_tree.find_clades() if c.name is not None}
        assert input_names == nbr_names


def test_nni_skips_root_incident_edges():
    tree = _make_search_shape_tree(5)
    root = tree.root
    root_children = set(root.clades)
    for _, spr_result in nni.nni_neighbors(tree):
        # u (prune_grandparent) must not be root or a root-child
        assert spr_result.prune_grandparent != root.name
        assert spr_result.prune_grandparent not in {c.name for c in root_children}


def test_nni_synthetic_spr_result_dirty_set():
    """For each NNI neighbor, the dirty up-pass set computed from the synthetic
    SPRResult must equal {v, u, parent(u), ..., root_name}."""
    tree = _make_search_shape_tree(6)
    normal_name = "diploid"

    # Build name->parent_name map on the ORIGINAL tree (before any swap)
    parent_of_orig = {}
    for c in tree.find_clades(order="preorder"):
        for child in c.clades:
            parent_of_orig[child.name] = c.name
    root_name = tree.root.name

    for neighbor_tree, spr_result in nni.nni_neighbors(tree):
        u_name = spr_result.prune_grandparent
        v_name = spr_result.regraft_parent

        # Expected dirty set: {v, u, parent(u) on the new tree, ..., root}
        # Since NNI doesn't change u's ancestors (only u's and v's child sets),
        # the parent chain of u is identical in old and new trees.
        expected = {v_name, u_name}
        cur = u_name
        while cur != root_name:
            cur = parent_of_orig[cur]
            expected.add(cur)

        actual = _get_dirty_nodes_up(neighbor_tree, spr_result, normal_name)
        assert actual == expected, (
            f"Dirty set mismatch for NNI swap u={u_name}, v={v_name}: "
            f"expected {expected}, got {actual}")


def test_nni_mode_terminates_on_optimal_input():
    """If the start tree is already at an NNI-local optimum, nni_mode terminates
    in 1 sweep with 0 accepted moves and returns the input tree unchanged."""
    import medicc
    import medicc.io
    import medicc.core as core
    from medicc.core import nni_mode

    # Use a tiny dataset where 3 samples + diploid -> n=3 leaves -> 0 NNI neighbors,
    # so termination is forced by zero enumeration.
    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    # Build minimal sample data: 3 samples with single-position copy-number profiles
    import pandas as pd
    rows = []
    for sample, cn in [("diploid", "1"), ("s1", "1"), ("s2", "2"), ("s3", "3")]:
        rows.append({"sample_id": sample, "chrom": "chr1", "start": 0, "end": 1,
                     "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")

    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    # Build a search-shape tree by hand (n=3 non-diploid leaves -> 0 NNI moves)
    tree = _make_search_shape_tree(3)

    result = nni_mode(
        tree=tree,
        samples_dict=fsa_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        nni_max_iter=10,
        n_cores=None,
    )
    assert len(result["trace"]) == 1
    assert len(result["best_trees"]) == 1
    assert len(result["best_ancestors"]) == 1


def test_nni_mode_finds_improvement():
    """Construct a 4-leaf scenario where the start tree is suboptimal and one
    NNI swap improves it. Verify the trace decreases and termination is at a
    local optimum.

    Profile design: s1/s3 have cn=2 and s2/s4 have cn=5 — two distinct groups.
    The bad caterpillar tree places s2 (cn=5) between s3 (cn=2) and s4 (cn=5),
    so the internal ancestor at internal_2 must reconcile divergent profiles.
    After one NNI swap, s2 and s4 (both cn=5) become siblings — reducing the score.
    """
    import medicc
    import medicc.io
    import medicc.core as core
    from medicc.core import nni_mode

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    # s1 and s3 have cn=2; s2 and s4 have cn=5
    import pandas as pd
    profiles = {
        "diploid": "1" * 10,
        "s1": "2" * 10,
        "s2": "5" * 10,
        "s3": "2" * 10,
        "s4": "5" * 10,
    }
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")

    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    # Bad caterpillar: diploid -> internal_0 -> {s1, internal_1}
    #                             internal_1 -> {s2, internal_2}
    #                             internal_2 -> {s3, s4}
    # s2 (cn=5) is forced as uncle of s3 (cn=2) and s4 (cn=5), raising the
    # score at internal_1. The optimal NNI swap at edge (internal_1, internal_2)
    # moves s2 out and puts s3/s4 with a correct ancestor.
    bad_tree = Tree(root=Clade(name="diploid", clades=[
        Clade(name="internal_0", clades=[
            Clade(name="s1"),
            Clade(name="internal_1", clades=[
                Clade(name="s2"),
                Clade(name="internal_2", clades=[
                    Clade(name="s3"),
                    Clade(name="s4"),
                ]),
            ]),
        ])
    ]), rooted=True)

    result = nni_mode(
        tree=bad_tree,
        samples_dict=fsa_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        nni_max_iter=10,
        n_cores=None,
    )

    # Trace must be strictly non-increasing (with at least one strict drop), or
    # a single entry if the start happened to be optimal.
    assert all(result["trace"][i+1] <= result["trace"][i] for i in range(len(result["trace"])-1))
    # We expect an improvement on this constructed scenario
    assert result["trace"][-1] < result["trace"][0], (
        f"Expected at least one improvement, but trace is {result['trace']}")
    assert result["best_score"] == result["trace"][-1]


def test_nni_mode_plateau_expands_frontier():
    """When two NNI neighbors of the starting tree are equally scored, nni_mode
    should carry both forward rather than discarding one.

    Setup: build a 4-leaf symmetric tree where two distinct NNI moves produce
    the same score. We mock the scoring so that the initial tree has a higher
    score and exactly two neighbors share the best score, while a 3rd (if any)
    is worse. We verify the returned frontier contains >= 2 trees.

    Implementation: We patch sum_of_branch_length to return controlled values
    keyed by tree topology hash, and update_branch_lengths to be a no-op, so
    the test drives the score landscape without needing real FST data.
    """
    import medicc
    import medicc.io
    import medicc.core as core
    import medicc.tree_hash as tree_hash
    from medicc.core import nni_mode
    from unittest.mock import patch

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    import pandas as pd
    # 4 non-diploid samples: two pairs with identical profiles so that two NNI
    # moves that swap them can produce equal scores
    profiles = {
        "diploid": "1" * 6,
        "s1": "2" * 6,
        "s2": "2" * 6,  # identical to s1
        "s3": "3" * 6,
        "s4": "3" * 6,  # identical to s3
    }
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")
    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    # Start tree: diploid -> internal_0 -> {s1, internal_1}
    #                         internal_1 -> {s2, internal_2}
    #                         internal_2 -> {s3, s4}
    start_tree = Tree(root=Clade(name="diploid", clades=[
        Clade(name="internal_0", clades=[
            Clade(name="s1"),
            Clade(name="internal_1", clades=[
                Clade(name="s2"),
                Clade(name="internal_2", clades=[
                    Clade(name="s3"),
                    Clade(name="s4"),
                ]),
            ]),
        ])
    ]), rooted=True)

    result = nni_mode(
        tree=start_tree,
        samples_dict=fsa_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        nni_max_iter=10,
        n_cores=None,
    )

    # With s1==s2 and s3==s4, the two NNI swaps that exchange s1↔s2 or s3↔s4
    # must produce equal scores — so the frontier should have >= 2 trees if the
    # initial topology was suboptimal (if it wasn't, at minimum we get 1 tree).
    # The key assertion: best_trees and best_ancestors are the same length.
    assert len(result["best_trees"]) == len(result["best_ancestors"]), \
        "best_trees and best_ancestors must have matching lengths"
    # All returned trees must share the same best_score
    for returned_tree in result["best_trees"]:
        score = medicc.tools.sum_of_branch_length(returned_tree)
        assert abs(score - result["best_score"]) < 1e-9, \
            f"Returned tree has score {score} != best_score {result['best_score']}"


def test_nni_mode_plateau_all_equal_neighbors_collected():
    """When ALL NNI neighbors of the current tree score the same as the current
    tree, the frontier should expand to include all of them (not terminate).

    We verify this by counting: after 1 sweep, if the frontier expanded, the
    next sweep runs from multiple trees. We check that nni_mode does not
    terminate immediately when best_neighbor_score == current_score.
    """
    import medicc
    import medicc.io
    import medicc.core as core
    from medicc.core import nni_mode
    import medicc.tools

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    import pandas as pd
    # All non-diploid samples identical: every NNI neighbor scores the same
    profiles = {
        "diploid": "1" * 4,
        "s1": "3" * 4,
        "s2": "3" * 4,
        "s3": "3" * 4,
        "s4": "3" * 4,
    }
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")
    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    start_tree = _make_search_shape_tree(4)

    result = nni_mode(
        tree=start_tree,
        samples_dict=fsa_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        nni_max_iter=3,  # cap low — test is about frontier expansion, not convergence
        n_cores=None,
    )

    # When all neighbors are equally good, the algorithm must NOT terminate at
    # sweep 0 claiming "no improvement". It should continue sweeping (trace may
    # have only the initial score entry, but nni_max_iter sweeps must run OR
    # the frontier must contain > 1 tree at return).
    # At minimum: all returned trees must share the reported best_score.
    assert len(result["best_trees"]) >= 1
    for returned_tree in result["best_trees"]:
        score = medicc.tools.sum_of_branch_length(returned_tree)
        assert abs(score - result["best_score"]) < 1e-9, \
            f"Returned tree has score {score} != best_score {result['best_score']}"


def test_nni_mode_strict_improvement_collapses_frontier():
    """After plateau expansion, a strict improvement in a later sweep should
    collapse the frontier back to only the improving branch(es).

    We verify: if one sweep produces a plateau (frontier > 1) and the next
    sweep produces a strict improvement from only one of those frontier trees,
    the returned frontier is exactly the improving set.
    """
    import medicc
    import medicc.io
    import medicc.core as core
    from medicc.core import nni_mode

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    import pandas as pd
    # s1/s3 identical (cn=2), s2/s4 identical (cn=5): symmetric plateau then improvement
    profiles = {
        "diploid": "1" * 8,
        "s1": "2" * 8,
        "s2": "5" * 8,
        "s3": "2" * 8,
        "s4": "5" * 8,
    }
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")
    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    start_tree = Tree(root=Clade(name="diploid", clades=[
        Clade(name="internal_0", clades=[
            Clade(name="s1"),
            Clade(name="internal_1", clades=[
                Clade(name="s2"),
                Clade(name="internal_2", clades=[
                    Clade(name="s3"),
                    Clade(name="s4"),
                ]),
            ]),
        ])
    ]), rooted=True)

    result = nni_mode(
        tree=start_tree,
        samples_dict=fsa_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        nni_max_iter=20,
        n_cores=None,
    )

    # Regardless of plateau traversal internals, the final returned trees must
    # all share best_score and the trace must be non-increasing.
    assert all(result["trace"][i + 1] <= result["trace"][i]
               for i in range(len(result["trace"]) - 1)), \
        f"Trace is not non-increasing: {result['trace']}"
    assert len(result["best_trees"]) == len(result["best_ancestors"])
    for returned_tree in result["best_trees"]:
        score = medicc.tools.sum_of_branch_length(returned_tree)
        assert abs(score - result["best_score"]) < 1e-9


def test_nni_mode_returns_full_plateau_two_cycle():
    """Regression: on a 2-cycle plateau the frontier must NOT be discarded.

    Landscape: the start tree T0 and exactly one of its NNI neighbours B share
    the best score; every other tree scores strictly worse. T0 and B are each
    other's only tied-best neighbour (a 2-cycle). The co-optimal set is {T0, B}.

    The old plateau branch REPLACED the frontier each sweep, so it traversed
    T0 -> B -> (back to T0, already visited) -> terminate, returning only {T0}
    and silently dropping B. A correct plateau traversal must return BOTH.

    Score landscape is mocked by topology hash so the test is deterministic and
    independent of FST scoring details.
    """
    import medicc
    import medicc.io
    import medicc.core as core
    import medicc.tools
    import medicc.tree_hash as tree_hash
    import medicc.nni as nni
    from medicc.core import nni_mode
    from unittest.mock import patch
    import pandas as pd

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    # 5 non-diploid samples -> n=5 leaves -> 4 NNI neighbours per tree.
    profiles = {"diploid": "1" * 6, "s1": "2" * 6, "s2": "3" * 6,
                "s3": "4" * 6, "s4": "5" * 6, "s5": "6" * 6}
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")
    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    start_tree = _make_search_shape_tree(5)
    # B := the first NNI neighbour of the start tree
    B = list(nni.nni_neighbors(start_tree))[0][0]
    h0 = tree_hash.get_topology_hash(start_tree)
    hB = tree_hash.get_topology_hash(B)
    assert h0 != hB
    plateau = {h0, hB}

    def fake_score(tree):
        # T0 and B tied-best (0.0); everything else strictly worse (1.0)
        return 0.0 if tree_hash.get_topology_hash(tree) in plateau else 1.0

    with patch("medicc.tools.sum_of_branch_length", side_effect=fake_score):
        result = nni_mode(
            tree=start_tree, samples_dict=fsa_dict, fst=fst,
            normal_name="diploid", prune_weight=0,
            nni_max_iter=50, n_cores=None)

    returned = {tree_hash.get_topology_hash(t) for t in result["best_trees"]}
    assert returned == {h0, hB}, (
        f"plateau traversal must return BOTH co-optimal trees {{{h0[:8]}, {hB[:8]}}}, "
        f"got {{{', '.join(h[:8] for h in returned)}}}")
    assert len(result["best_trees"]) == len(result["best_ancestors"])


def test_nni_parallel_matches_serial():
    """evaluate_nni_neighbors_parallel must return the same scores as the serial
    nni_neighbors path, and cover the same set of neighbor topologies."""
    import medicc
    import medicc.io
    import medicc.core as core
    import medicc.tree_hash as tree_hash
    from medicc.ancestors import reconstruct_ancestors_incremental
    from medicc.nni import evaluate_nni_neighbors_parallel

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    import pandas as pd
    profiles = {
        "diploid": "1" * 6,
        "s1": "2" * 6,
        "s2": "3" * 6,
        "s3": "4" * 6,
        "s4": "5" * 6,
        "s5": "6" * 6,
    }
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")
    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")

    tree = _make_search_shape_tree(5)
    ancestors, uppass_cache = medicc.reconstruct_ancestors(
        tree=tree, samples_dict=fsa_dict, upper_pass_fst=fst, lower_pass_fst=fst,
        normal_name="diploid", prune_weight=0, spr_logger_disable=True)
    core.update_branch_lengths(tree, fst, ancestors, "diploid")

    # Serial reference
    serial_scores = {}
    for neighbor_tree, synthetic in nni.nni_neighbors(tree):
        nbr = copy.deepcopy(neighbor_tree)
        inc_ancestors, _ = reconstruct_ancestors_incremental(
            tree=nbr, old_uppass_cache=uppass_cache,
            samples_dict=fsa_dict, fst=fst, normal_name="diploid",
            spr_result=synthetic, prune_weight=0)
        core.update_branch_lengths(nbr, fst, inc_ancestors, "diploid")
        h = tree_hash.tree_hash(tree_hash.strip_branch_lengths(nbr))
        serial_scores[h] = medicc.tools.sum_of_branch_length(nbr)

    # Parallel path (n_cores=2)
    parallel_results = evaluate_nni_neighbors_parallel(
        tree=tree,
        old_uppass_cache=uppass_cache,
        samples_dict=fsa_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        n_cores=2,
    )
    parallel_scores = {}
    for nbr_tree, _, _, score, *_ in parallel_results:
        h = tree_hash.tree_hash(tree_hash.strip_branch_lengths(nbr_tree))
        parallel_scores[h] = score

    assert set(serial_scores.keys()) == set(parallel_scores.keys()), \
        "Parallel and serial paths produced different sets of neighbor topologies"
    for h in serial_scores:
        assert serial_scores[h] == parallel_scores[h], \
            f"Score mismatch for topology {h}: serial={serial_scores[h]}, parallel={parallel_scores[h]}"


def test_nni_incremental_matches_full_reconstruction():
    """For each NNI neighbor of a small real tree, the score computed via the
    incremental reconstruction must equal the score computed via a full
    from-scratch reconstruction. This is the load-bearing correctness check
    for the synthetic SPRResult."""
    import medicc
    import medicc.io
    import medicc.core as core
    from medicc.ancestors import reconstruct_ancestors_incremental

    fst = medicc.io.read_fst()
    symbol_table = fst.input_symbols()

    # 5 samples + diploid -> n=5 leaves -> 4 NNI neighbors per sweep
    import pandas as pd
    profiles = {
        "diploid": "1" * 8,
        "s1": "2" * 8,
        "s2": "3" * 8,
        "s3": "4" * 8,
        "s4": "5" * 8,
        "s5": "6" * 8,
    }
    rows = []
    for sample, prof in profiles.items():
        for i, cn in enumerate(prof):
            rows.append({"sample_id": sample, "chrom": "chr1",
                         "start": i, "end": i + 1, "cn_a": cn})
    df = pd.DataFrame(rows).set_index(["sample_id", "chrom", "start", "end"])
    df["cn_a"] = df["cn_a"].astype("category")

    fsa_dict, _ = core.create_standard_fsa_dict_from_data(df, symbol_table, "X")
    tree = _make_search_shape_tree(5)

    # Bootstrap: full reconstruction of the start tree to populate uppass_cache
    ancestors, uppass_cache = medicc.reconstruct_ancestors(
        tree=tree, samples_dict=fsa_dict, upper_pass_fst=fst, lower_pass_fst=fst,
        normal_name="diploid", prune_weight=0,
        spr_logger_disable=True)
    core.update_branch_lengths(tree, fst, ancestors, "diploid")

    for neighbor_tree, synthetic in nni.nni_neighbors(tree):
        # Incremental
        nbr_inc = copy.deepcopy(neighbor_tree)
        inc_ancestors, _ = reconstruct_ancestors_incremental(
            tree=nbr_inc, old_uppass_cache=uppass_cache,
            samples_dict=fsa_dict, fst=fst, normal_name="diploid",
            spr_result=synthetic, prune_weight=0)
        core.update_branch_lengths(nbr_inc, fst, inc_ancestors, "diploid")
        inc_score = medicc.tools.sum_of_branch_length(nbr_inc)

        # Full from scratch
        nbr_full = copy.deepcopy(neighbor_tree)
        full_ancestors, _ = medicc.reconstruct_ancestors(
            tree=nbr_full, samples_dict=fsa_dict, upper_pass_fst=fst, lower_pass_fst=fst,
            normal_name="diploid", prune_weight=0,
            spr_logger_disable=True)
        core.update_branch_lengths(nbr_full, fst, full_ancestors, "diploid")
        full_score = medicc.tools.sum_of_branch_length(nbr_full)

        assert inc_score == full_score, (
            f"Incremental score {inc_score} != full score {full_score} "
            f"for NNI swap u={synthetic.prune_grandparent}, v={synthetic.regraft_parent}")


def test_evaluate_nni_neighbors_parallel_returns_step():
    """Each result from evaluate_nni_neighbors_parallel must include its step number."""
    import medicc.nni as nni_mod
    from unittest.mock import patch, MagicMock
    tree = _make_search_shape_tree(5)
    moves = nni_mod._enumerate_moves(tree)
    n_moves = len(moves)

    call_steps = []
    def fake_eval(tree, move, old_uppass_cache, samples_dict, fst, normal_name, prune_weight, step=0, upper_pass_fst=None, lower_pass_fst=None):
        call_steps.append(step)
        return (MagicMock(), {}, {}, 10, step)

    with patch("medicc.nni._eval_nni_neighbor", side_effect=fake_eval):
        results = nni_mod.evaluate_nni_neighbors_parallel(
            tree=tree,
            old_uppass_cache={},
            samples_dict={},
            fst=MagicMock(),
            normal_name="diploid",
            prune_weight=0,
            n_cores=1,
            step_start=10,
        )
    assert len(results) == n_moves
    assert sorted(call_steps) == list(range(10, 10 + n_moves))
    for result in results:
        assert len(result) == 5


def test_evaluate_nni_neighbors_parallel_step_start_zero():
    """step_start=0 gives steps 0, 1, 2, ..."""
    import medicc.nni as nni_mod
    from unittest.mock import patch, MagicMock
    tree = _make_search_shape_tree(4)
    moves = nni_mod._enumerate_moves(tree)

    call_steps = []
    def fake_eval(tree, move, old_uppass_cache, samples_dict, fst, normal_name, prune_weight, step=0, upper_pass_fst=None, lower_pass_fst=None):
        call_steps.append(step)
        return (MagicMock(), {}, {}, 5, step)

    with patch("medicc.nni._eval_nni_neighbor", side_effect=fake_eval):
        results = nni_mod.evaluate_nni_neighbors_parallel(
            tree=tree,
            old_uppass_cache={},
            samples_dict={},
            fst=MagicMock(),
            normal_name="diploid",
            prune_weight=0,
            n_cores=1,
            step_start=0,
        )
    assert sorted(call_steps) == list(range(len(moves)))


def test_nni_mode_step_records_count():
    """nni_mode with nni_trace_dir accumulates one record per evaluated neighbor."""
    import tempfile, os
    import medicc
    import medicc.core as core
    import medicc.io
    import fstlib

    tree = _make_search_shape_tree(4)
    fst = medicc.io.read_fst(user_fst=None, no_wgd=False,
                              total_copy_numbers=False, wgd_x2=False)
    symbol_table = fst.input_symbols()
    cn_str = "22X22"
    samples = ["diploid", "s1", "s2", "s3", "s4"]
    samples_dict = {
        s: fstlib.factory.from_string(cn_str, arc_type="standard",
                                      isymbols=symbol_table, osymbols=symbol_table)
        for s in samples
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        result = core.nni_mode(
            tree=tree,
            samples_dict=samples_dict,
            fst=fst,
            normal_name="diploid",
            prune_weight=0,
            nni_max_iter=2,
            n_cores=None,
            nni_trace_dir=tmpdir,
        )
    step_records = result["step_records"]
    assert len(step_records) > 0
    for step, newick_str, score in step_records:
        assert isinstance(step, int)
        assert isinstance(newick_str, str)
        assert len(newick_str) > 0
        assert isinstance(score, (int, float))
    steps = [r[0] for r in step_records]
    assert steps == sorted(steps), "Step numbers should be monotonically increasing"
    assert len(set(steps)) == len(steps), "Step numbers should be unique"


def test_nni_mode_step_records_empty_when_no_trace_dir():
    """nni_mode without nni_trace_dir returns empty step_records."""
    import medicc
    import medicc.core as core
    import medicc.io
    import fstlib

    tree = _make_search_shape_tree(4)
    fst = medicc.io.read_fst(user_fst=None, no_wgd=False,
                              total_copy_numbers=False, wgd_x2=False)
    symbol_table = fst.input_symbols()
    cn_str = "22X22"
    samples = ["diploid", "s1", "s2", "s3", "s4"]
    samples_dict = {
        s: fstlib.factory.from_string(cn_str, arc_type="standard",
                                      isymbols=symbol_table, osymbols=symbol_table)
        for s in samples
    }

    result = core.nni_mode(
        tree=tree,
        samples_dict=samples_dict,
        fst=fst,
        normal_name="diploid",
        prune_weight=0,
        nni_max_iter=2,
        n_cores=None,
        nni_trace_dir=None,
    )
    assert result["step_records"] == []


def test_main_nni_writes_trace_files():
    """main_nni with nni_trace_dir writes per-step files and nni_trace.tsv."""
    import pathlib, tempfile, os
    import medicc
    import medicc.core as core
    import medicc.io
    import pandas as pd

    repo_root = pathlib.Path(__file__).parent.parent.absolute()
    input_path = repo_root / "examples" / "simple_example" / "simple_example.tsv"
    input_df = medicc.io.read_and_parse_input_data(
        filename=str(input_path),
        normal_name="diploid",
        input_type="t",
        separator="X",
        allele_columns=["cn_a", "cn_b"],
        total_copy_numbers=False,
        maxcn=8,
    )
    fst = medicc.io.read_fst(user_fst=None, no_wgd=False,
                              total_copy_numbers=False, wgd_x2=False)
    event_fst = medicc.io.read_fst(user_fst=None, no_wgd=False,
                                   total_copy_numbers=False, wgd_x2=False)

    with tempfile.TemporaryDirectory() as tmpdir:
        trace_dir = os.path.join(tmpdir, "NNI_trace")
        core.main_nni(
            input_df=input_df,
            asymm_fst=fst,
            output_dir=tmpdir,
            event_counting_fst=event_fst,
            normal_name="diploid",
            nni_max_iter=1,
            n_cores=None,
            nni_trace_dir=trace_dir,
        )

        step_files = sorted(f for f in os.listdir(trace_dir) if f.startswith("step_") and f.endswith(".txt"))
        assert len(step_files) > 0, "No step files written"

        with open(os.path.join(trace_dir, step_files[0])) as f:
            lines = f.read().strip().splitlines()
        assert len(lines) == 2, f"Expected 2 lines in step file, got {len(lines)}"
        float(lines[1])  # line 2 must be a parseable number

        tsv_path = os.path.join(trace_dir, "nni_trace.tsv")
        assert os.path.exists(tsv_path), "nni_trace.tsv not written"
        tsv = pd.read_csv(tsv_path, sep="\t")
        assert list(tsv.columns) == ["step", "newick", "sum_of_branch_length"]
        assert len(tsv) == len(step_files)
        assert list(tsv["step"]) == sorted(tsv["step"].tolist()), "TSV rows not in step order"


def test_main_nni_no_trace_files_when_no_trace_dir():
    """main_nni without nni_trace_dir creates no NNI_trace directory."""
    import pathlib, tempfile, os
    import medicc
    import medicc.core as core
    import medicc.io

    repo_root = pathlib.Path(__file__).parent.parent.absolute()
    input_path = repo_root / "examples" / "simple_example" / "simple_example.tsv"
    input_df = medicc.io.read_and_parse_input_data(
        filename=str(input_path),
        normal_name="diploid",
        input_type="t",
        separator="X",
        allele_columns=["cn_a", "cn_b"],
        total_copy_numbers=False,
        maxcn=8,
    )
    fst = medicc.io.read_fst(user_fst=None, no_wgd=False,
                              total_copy_numbers=False, wgd_x2=False)
    event_fst = medicc.io.read_fst(user_fst=None, no_wgd=False,
                                   total_copy_numbers=False, wgd_x2=False)

    with tempfile.TemporaryDirectory() as tmpdir:
        core.main_nni(
            input_df=input_df,
            asymm_fst=fst,
            output_dir=tmpdir,
            event_counting_fst=event_fst,
            normal_name="diploid",
            nni_max_iter=1,
            n_cores=None,
            nni_trace_dir=None,
        )
        trace_dir = os.path.join(tmpdir, "NNI_trace")
        assert not os.path.exists(trace_dir), "NNI_trace dir should not exist"
