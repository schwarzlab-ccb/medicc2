import logging
from collections import defaultdict
import copy 

import Bio
import Bio.Phylo
import fstlib
import numpy as np

import medicc

logger = logging.getLogger(__name__)


def _intersect_task(left_fsa, right_fsa, fst, prune_weight):
    """Worker function for parallel intersection (up-the-tree phase)."""
    return intersect_clades_detmin(
        left_fsa, right_fsa, fst,
        prune_weight=prune_weight,
        detmin_before_intersect=False,
        detmin_after_intersect=True,
    )


def _align_task(parent_fsa, child_fsa, fst):
    """Worker function for parallel alignment (down-the-tree phase)."""
    sp = fstlib.align(fst, parent_fsa, child_fsa)
    return fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')


def _group_nodes_by_depth(clade_list, normal_name):
    """Group internal nodes by depth in the tree.

    Given a preorder clade_list (excluding normal), assigns each node a depth
    from the root and returns a dict mapping depth -> list of nodes at that depth.
    Only internal nodes (nodes with children) are included.

    Note: Parallelism is limited by tree shape — a highly unbalanced (caterpillar)
    tree has at most 1 node per depth level, yielding no parallelism. Balanced
    binary trees benefit most.
    """
    depth_map = {id(clade_list[0]): 0}  # root is depth 0
    levels = defaultdict(list)

    for node in clade_list:
        node_depth = depth_map.get(id(node), 0)
        if len(node.clades) != 0:
            levels[node_depth].append(node)
            for child in node.clades:
                if child.name != normal_name:
                    depth_map[id(child)] = node_depth + 1

    return levels


def reconstruct_ancestors(tree, samples_dict, upper_pass_fst, normal_name, lower_pass_fst, prune_weight=0,
                          spr_logger_disable=False, n_cores=None):
    """Full ancestor reconstruction.

    upper_pass_fst: the FST to use during upper pass
    lower_pass_fst: the FST to use during lower pass

    Returns:
        tuple: (fsa_dict, uppass_cache)
            - fsa_dict: dict mapping node names to final (post-down-pass) FSAs
            - uppass_cache: dict mapping internal node names to their up-pass
              (intersection) FSAs, before the down-pass overwrites them.
              This is used by reconstruct_ancestors_incremental to cache
              clean nodes across SPR iterations.
    """
    # spr_logger: disable logging for ancestor reconstruction in spr mode for simplicity
    if len(samples_dict) == 2:
        return samples_dict, {}

    fsa_dict = samples_dict.copy()
    tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)

    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]

    original_logger_level = logger.getEffectiveLevel()
    original_propagate = logger.propagate
    original_handlers = logger.handlers.copy()

    if spr_logger_disable:
        logger.setLevel(logging.CRITICAL + 1)
        logger.propagate = False  # <-- Prevent bubbling up
        logger.handlers = [logging.NullHandler()]

    use_parallel = n_cores is not None and n_cores > 1
    levels = _group_nodes_by_depth(clade_list, normal_name)

    logger.info("Ancestor reconstruction: Up the tree")
    # up the tree (leaf to root), level by level (bottom-up)
    for depth in sorted(levels.keys(), reverse=True):
        nodes_at_depth = [n for n in levels[depth] if len(n.clades) != 0]
        if not nodes_at_depth:
            continue

        if use_parallel and len(nodes_at_depth) > 1:
            from joblib import Parallel, delayed
            results = Parallel(n_jobs=n_cores)(
                delayed(_intersect_task)(
                    fsa_dict[[c for c in node.clades if c.name != normal_name][0].name],
                    fsa_dict[[c for c in node.clades if c.name != normal_name][1].name],
                    upper_pass_fst, prune_weight
                ) for node in nodes_at_depth
            )
            for node, result in zip(nodes_at_depth, results):
                logger.debug(f"Clade: {node.name} (parallel)")
                fsa_dict[node.name] = result
        else:
            for node in nodes_at_depth:
                children = [item for item in node.clades if item.name != normal_name]
                left_name = children[0].name
                right_name = children[1].name
                logger.debug(f"Clade: {node.name}, left: {left_name}, right: {right_name}")
                fsa_dict[node.name] = _intersect_task(
                    fsa_dict[left_name], fsa_dict[right_name], upper_pass_fst, prune_weight)

    # Snapshot up-pass results before root alignment and down-pass overwrite them
    uppass_cache = {node.name: fsa_dict[node.name]
                    for node in clade_list if len(node.clades) != 0}

    logger.debug("Ancestor reconstruction for root")
    # root node is calculated separately w.r.t. normal node
    root_name = clade_list[0].name 
    sp = fstlib.align(lower_pass_fst, fsa_dict[normal_name], fsa_dict[root_name])
    fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    logger.info("Ancestor reconstruction: Down the tree")
    # down the tree (root to leaf), level by level (top-down)
    for depth in sorted(levels.keys()):
        pairs = []
        for node in levels[depth]:
            if len(node.clades) != 0:
                for child in node.clades:
                    if len(child.clades) != 0:
                        pairs.append((node.name, child.name))

        if not pairs:
            continue

        if use_parallel and len(pairs) > 1:
            from joblib import Parallel, delayed
            results = Parallel(n_jobs=n_cores)(
                delayed(_align_task)(fsa_dict[pname], fsa_dict[cname], lower_pass_fst)
                for pname, cname in pairs
            )
            for (_, cname), result in zip(pairs, results):
                logger.debug(f"Aligned child: {cname} (parallel)")
                fsa_dict[cname] = result
        else:
            for pname, cname in pairs:
                logger.debug(f"Clade: {pname}, internal child: {cname}")
                fsa_dict[cname] = _align_task(fsa_dict[pname], fsa_dict[cname], lower_pass_fst)

    # check if ancestors were correctly reconstructed
    sample_lengths = {sample: len(medicc.tools.fsa_to_string(fsa_dict[sample])) for sample, fsa in fsa_dict.items()}
    normal_length = sample_lengths[normal_name]

    if np.any([x != normal_length for x in sample_lengths.values()]):
        raise MEDICCAncestorReconstructionError("Some ancestors could not be reconstructed. These are:\n"
                                                "{}".format('\n'.join([sample for sample, length in sample_lengths.items() if length != normal_length])) + \
                                                "\nCheck whether your normal sample contains segments with copy number zero")

    if spr_logger_disable:
        logger.setLevel(original_logger_level)
        logger.propagate = original_propagate
        logger.handlers = original_handlers
    return fsa_dict, uppass_cache


def _detmin_weighted(f):
    """Determinize+minimize a weighted acceptor to keep the cost-to-go automata
    compact (the project step below creates one parallel path per child state;
    determinization collapses them to one min-weight path per string). Weighted
    determinization is exact for these acyclic profile automata."""
    return fstlib.determinize(f.arcsort('ilabel')).minimize()


def _cost_to_go_through_child(obj_fst, child_cost_acc, parent_mask):
    """Cost-to-go contribution of one child, as a function of the parent state s:
        T(s) = min over child states s_c of [ obj_fst(s -> s_c) + C_child(s_c) ]
    masked to the parent's event-minimal candidate set. The project to the input
    side performs the min over s_c in the tropical semiring."""
    comp = fstlib.compose(obj_fst.arcsort('olabel'), child_cost_acc.copy().arcsort('ilabel'))
    T = comp.project('input')
    T = fstlib.intersect(T.arcsort('olabel'), parent_mask)
    return _detmin_weighted(T)


def _best_child(obj_fst, parent_string_fsa, child_cost_acc):
    """Cost-aware backtrace for one parent->child edge: given the already-committed
    single parent string, return the child string minimising
        obj_fst(parent -> s_c) + C_child(s_c)."""
    comp = fstlib.compose(parent_string_fsa.copy().arcsort('olabel'), obj_fst.arcsort('olabel'))
    comp = fstlib.compose(comp.arcsort('olabel'), child_cost_acc.copy().arcsort('ilabel'))
    sp = fstlib.shortestpath(comp)
    return fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')


def reconstruct_ancestors_constrained_sankoff(tree, samples_dict, upper_pass_fst, normal_name,
                                              lower_pass_fst, prune_weight=0,
                                              spr_logger_disable=False, n_cores=None):
    """Constrained-Sankoff ancestor reconstruction.

    Same up-pass as reconstruct_ancestors (the ``upper_pass_fst`` builds each
    internal node's *event-minimal* candidate set), but the greedy parent-only
    down-pass is replaced by a proper Sankoff dynamic program under
    ``lower_pass_fst``: a cost-to-go up-sweep (best achievable objective cost over
    the whole subtree as a function of the node's state, restricted to the
    event-minimal candidate set) followed by a cost-aware backtrace.

    With ``lower_pass_fst`` an open-K length-encoding FST (K larger than the
    maximum total span, e.g. open=5000), this returns the reconstruction that
    minimises the total focal-event length among ALL event-minimal labellings
    (whole-chromosome / WGD events are flat-cost and contribute no length). Unlike
    the greedy down-pass it accounts for each node's descendants, so it never does
    worse than greedy and is provably optimal for the objective over the
    event-minimal feasible set.

    Args mirror reconstruct_ancestors; returns (fsa_dict, uppass_cache).

    Note: the cost-to-go automata can be large for trees with very many co-optimal
    ancestors; this runs serially and is memory-bound on such cases (n_cores is
    accepted for signature compatibility but unused).
    """
    if len(samples_dict) == 2:
        return samples_dict, {}

    fsa_dict = samples_dict.copy()
    tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)
    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]
    levels = _group_nodes_by_depth(clade_list, normal_name)

    # --- up the tree: event-minimal candidate sets S_v (event-counting FST) ---
    logger.info("Constrained Sankoff: up the tree (event-minimal candidate sets)")
    for depth in sorted(levels.keys(), reverse=True):
        for node in (n for n in levels[depth] if len(n.clades) != 0):
            children = [item for item in node.clades if item.name != normal_name]
            fsa_dict[node.name] = intersect_clades_detmin(
                fsa_dict[children[0].name], fsa_dict[children[1].name], upper_pass_fst,
                prune_weight=prune_weight, detmin_before_intersect=False, detmin_after_intersect=True)

    uppass_cache = {node.name: fsa_dict[node.name]
                    for node in clade_list if len(node.clades) != 0}

    # --- cost-to-go under the objective FST, masked to S_v at every node ---
    logger.info("Constrained Sankoff: cost-to-go up-sweep (objective FST)")
    cost_to_go = {name: samples_dict[name] for name in samples_dict}  # leaf cost-to-go = leaf string, weight 0
    for depth in sorted(levels.keys(), reverse=True):
        for node in (n for n in levels[depth] if len(n.clades) != 0):
            children = [item for item in node.clades if item.name != normal_name]
            mask = fstlib.arcmap(fsa_dict[node.name].copy(), map_type='rmweight').arcsort('ilabel')
            terms = [_cost_to_go_through_child(lower_pass_fst, cost_to_go[c.name], mask) for c in children]
            Cv = terms[0]
            for t in terms[1:]:
                Cv = _detmin_weighted(fstlib.intersect(Cv.arcsort('olabel'), t.arcsort('ilabel')))
            cost_to_go[node.name] = Cv

    # --- root choice, then cost-aware backtrace down the tree ---
    logger.info("Constrained Sankoff: backtrace down the tree")
    root_name = clade_list[0].name
    fsa_dict[root_name] = _best_child(lower_pass_fst, fsa_dict[normal_name], cost_to_go[root_name])
    for depth in sorted(levels.keys()):
        for node in levels[depth]:
            if len(node.clades) == 0:
                continue
            for child in node.clades:
                if child.name == normal_name or len(child.clades) == 0:
                    continue
                fsa_dict[child.name] = _best_child(lower_pass_fst, fsa_dict[node.name], cost_to_go[child.name])

    # check if ancestors were correctly reconstructed
    sample_lengths = {sample: len(medicc.tools.fsa_to_string(fsa_dict[sample])) for sample, fsa in fsa_dict.items()}
    normal_length = sample_lengths[normal_name]
    if np.any([x != normal_length for x in sample_lengths.values()]):
        raise MEDICCAncestorReconstructionError("Some ancestors could not be reconstructed. These are:\n"
                                                "{}".format('\n'.join([sample for sample, length in sample_lengths.items() if length != normal_length])) + \
                                                "\nCheck whether your normal sample contains segments with copy number zero")
    return fsa_dict, uppass_cache


def _get_dirty_nodes_up(tree, spr_result, normal_name):
    """Compute the set of internal node names whose up-pass (intersection) results
    need recomputation after an SPR move.

    The dirty set is the union of two root-ward paths:
      1. From the prune grandparent (G) up to the root.
      2. From the regraft parent (P') up to the root.

    These are the only nodes whose descendant leaf sets changed.
    """
    dirty = set()
    root_name = tree.root.name

    # Build a child->parent lookup (by name) for the NEW tree
    parent_of = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            if child.name != normal_name:
                parent_of[child.name] = clade.name

    # Walk from prune_grandparent to root
    node_name = spr_result.prune_grandparent
    while node_name is not None:
        dirty.add(node_name)
        if node_name == root_name:
            break
        node_name = parent_of.get(node_name)

    # Walk from regraft_parent to root
    node_name = spr_result.regraft_parent
    while node_name is not None:
        dirty.add(node_name)
        if node_name == root_name:
            break
        node_name = parent_of.get(node_name)

    return dirty


def reconstruct_ancestors_incremental(tree, old_uppass_cache, samples_dict, fst, normal_name,
                                      spr_result, prune_weight=0):
    """Incremental ancestor reconstruction after an SPR move.

    Reuses the up-pass (intersection) results from old_uppass_cache for nodes whose
    descendant leaf sets did not change. Only recomputes up-pass for 'dirty' nodes
    (those on the path from prune/regraft points to root).

    The down-pass (alignment against parent) is always run in full, since changes
    propagate from the root downward and it's simpler/safer to redo entirely.

    Args:
        old_uppass_cache: dict mapping internal node names to their up-pass
            (intersection-only) FSAs from the previous iteration. This must be
            the pure up-pass results, NOT the post-down-pass aligned ancestors.

    Returns:
        tuple: (fsa_dict, uppass_cache) — same contract as reconstruct_ancestors.
    """
    fsa_dict = samples_dict.copy()
    tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)

    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]

    if len(clade_list) == 0:
        return fsa_dict, {}

    # Determine which nodes need up-pass recomputation
    dirty_up = _get_dirty_nodes_up(tree, spr_result, normal_name)

    logger.debug(f"Incremental reconstruction: {len(dirty_up)} dirty up-pass nodes out of {len(clade_list)}")

    # --- UP THE TREE (leaf to root) ---
    # Reuse old up-pass results for clean nodes, recompute for dirty ones
    for node in reversed(clade_list):
        if len(node.clades) != 0:
            if node.name in dirty_up:
                # Recompute intersection for this dirty node
                children = [item for item in node.clades if item.name != normal_name]
                left_name = children[0].name
                right_name = children[1].name
                logger.debug(f"Incremental up-pass (dirty): {node.name}, left: {left_name}, right: {right_name}")
                fsa_dict[node.name] = _intersect_task(
                    fsa_dict[left_name], fsa_dict[right_name], fst, prune_weight)
            else:
                # Reuse cached up-pass intersection result
                if node.name in old_uppass_cache:
                    logger.debug(f"Incremental up-pass (cached): {node.name}")
                    fsa_dict[node.name] = old_uppass_cache[node.name]
                else:
                    # Node name exists in new tree but not in old cache (shouldn't happen
                    # since SPR reuses internal node names, but handle gracefully)
                    children = [item for item in node.clades if item.name != normal_name]
                    left_name = children[0].name
                    right_name = children[1].name
                    logger.debug(f"Incremental up-pass (new node): {node.name}, left: {left_name}, right: {right_name}")
                    fsa_dict[node.name] = _intersect_task(
                        fsa_dict[left_name], fsa_dict[right_name], fst, prune_weight)

    # Snapshot up-pass results before root alignment and down-pass overwrite them
    uppass_cache = {node.name: fsa_dict[node.name]
                    for node in clade_list if len(node.clades) != 0}

    # --- ROOT ---
    root_name = clade_list[0].name
    sp = fstlib.align(fst, fsa_dict[normal_name], fsa_dict[root_name])
    fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    # --- DOWN THE TREE (root to leaf) — always full ---
    for node in clade_list:
        if len(node.clades) != 0:
            children = [q for q in node.clades if len(q.clades) != 0]
            for child in children:
                sp = fstlib.align(fst, fsa_dict[node.name], fsa_dict[child.name])
                fsa_dict[child.name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    # Validation
    sample_lengths = {sample: len(medicc.tools.fsa_to_string(fsa_dict[sample])) for sample, fsa in fsa_dict.items()}
    normal_length = sample_lengths[normal_name]

    if np.any([x != normal_length for x in sample_lengths.values()]):
        raise MEDICCAncestorReconstructionError("Some ancestors could not be reconstructed. These are:\n"
                                                "{}".format('\n'.join([sample for sample, length in sample_lengths.items() if length != normal_length])) + \
                                                "\nCheck whether your normal sample contains segments with copy number zero")

    return fsa_dict, uppass_cache


def intersect_clades_detmin(left, right, fst, prune_weight=None, detmin_before_intersect=True, detmin_after_intersect=True):
    L = fstlib.compose(fst, left.arcsort('ilabel')).project('input')
    R = fstlib.compose(fst, right.arcsort('ilabel')).project('input')
    if detmin_before_intersect:
        L = fstlib.determinize(L).minimize()
        R = fstlib.determinize(R).minimize()
    intersection = fstlib.intersect(L.arcsort('olabel'), R)
    # For prune_weight=0, deletes all paths but the shortest one
    if prune_weight is not None:
        pruned = fstlib.prune(intersection, weight=prune_weight)
    else:
        pruned = intersection
    if detmin_after_intersect:
        pruned = fstlib.determinize(pruned).minimize()
    return pruned


class MEDICCAncestorReconstructionError(Exception):
    pass
