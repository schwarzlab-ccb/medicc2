import logging
from collections import defaultdict
from joblib import Parallel, delayed

import Bio
import Bio.Phylo
import fstlib
import numpy as np

import medicc

logger = logging.getLogger(__name__)

def _intersect_task(left_fsa, right_fsa, fst, prune_weight=None, detmin_before_intersect=False, detmin_after_intersect=True):
    """Worker function for parallel intersection (up-the-tree)"""
    return intersect_clades_detmin(
        left_fsa, right_fsa, fst,
        prune_weight=prune_weight,
        detmin_before_intersect=detmin_before_intersect,
        detmin_after_intersect=detmin_after_intersect,
    )

def _align_task(parent_fsa, child_fsa, fst):
    """Worker function for parallel alignment (down-the-tree)"""
    sp = fstlib.align(fst, parent_fsa, child_fsa)
    return fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

def _group_nodes_by_depth(clade_list, normal_name):
    """
    Group internal nodes by depth in the tree

    Given a preorder claide_list (excluding normal), assigns each node a depth from the root and returns a dict mapping
    depth -> list of nodes that depth. Only internal nodes are included.

    Used for parallel ancestral reconstruction, and parallelism is limited by tree shape -- a highly unbalanced (caterpillar)
    tree has at most 1 node per depth level, yielding no parallelism. Balanced binary trees benefit most.
    """

    depth_map = {id(clade_list[0]): 0}
    levels = defaultdict(list)

    # breadth first search of the tree.
    for node in clade_list:
        node_depth = depth_map.get(id(node), 0) # default 0
        if len(node.clades) != 0:
            levels[node_depth].append(node)
            for child in node.clades:
                if child.name != normal_name:
                    depth_map[id(child)] = node_depth + 1
    return levels



def reconstruct_ancestors(tree, samples_dict, upper_pass_fst, lower_pass_fst, normal_name, prune_weight=0, n_cores=None, upper_cache=False):
    '''
    Reconstruct ancestors: up-pass intersects sibling FSAs under upper_pass_fst,
    then a down-pass aligns root-to-normal and parent-to-child under lower_pass_fst.
    Does not mutate `tree` or `samples_dict`.

    n_cores: None or 1 runs sequentially; any other value parallelizes each
    tree-depth level with joblib. Benefit depends on tree shape — a caterpillar
    tree has ~1 node per depth and gains nothing.

    upper_cache: when True, also returns uppass_cache, the pre-down-pass
    up-pass (intersection-only) FSA for each internal node — kept for callers
    that later want to reuse unchanged subtrees instead of recomputing them.

    Returns fsa_dict, or (fsa_dict, uppass_cache) if upper_cache=True.
    '''

    if len(samples_dict) == 2:
        if not upper_cache:
            return samples_dict
        else:
            return samples_dict, {}

    fsa_dict = samples_dict.copy()
    tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)

    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]

    if n_cores is None or n_cores == 1:
        # Sequential pipeline
        logger.info("Ancestor reconstruction: Up the tree")
        # up the tree (leaf to root)
        for node in reversed(clade_list):
            if len(node.clades) != 0:
                children = [item for item in node.clades if item.name != normal_name]
                left_name = children[0].name
                right_name = children[1].name
                logger.debug(f"Clade: {node.name}, left: {left_name}, right: {right_name}")

                ## project
                intersection = intersect_clades_detmin(fsa_dict[left_name], fsa_dict[right_name], upper_pass_fst,
                                                       prune_weight=prune_weight, detmin_before_intersect=False, detmin_after_intersect=True)
                fsa_dict[node.name] = intersection

        if upper_cache:
            uppass_cache = {node.name: fsa_dict[node.name]
                            for node in clade_list if len(node.clades) != 0}

        logger.debug("Ancestor reconstruction for root")
        # root node is calculated separately w.r.t. normal node
        root_name = clade_list[0].name
        sp = fstlib.align(lower_pass_fst, fsa_dict[normal_name], fsa_dict[root_name])
        fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

        logger.info("Ancestor reconstruction: Down the tree")
        # down the tree (root to leaf)
        for node in clade_list:
            if len(node.clades) != 0:
                children = [q for q in node.clades if len(q.clades) != 0]
                logger.debug(f"Clade: {node.name}, internal children: {children}")
                for child in children:
                    sp = fstlib.align(lower_pass_fst, fsa_dict[node.name], fsa_dict[child.name])
                    fsa_dict[child.name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')
    else:
        # parallel pipeline
        logger.info(f"Ancestor reconstruction using {n_cores} threads: Up the tree")
        levels = _group_nodes_by_depth(clade_list, normal_name)
        for depth in sorted(levels.keys(), reverse=True):
            nodes_at_depth = [n for n in levels[depth] if len(n.clades) != 0]
            if not nodes_at_depth:
                continue


            if len(nodes_at_depth) > 1: # if one node you pay parallize overhead for nothing
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
                    fsa_dict[node.name] = _intersect_task(fsa_dict[left_name], fsa_dict[right_name], upper_pass_fst, prune_weight)

        if upper_cache:
            uppass_cache = {node.name: fsa_dict[node.name]
                            for node in clade_list if len(node.clades) != 0}

        logger.debug("Ancestor reconstruction for root")
        # root node is calculated separately w.r.t. normal node
        root_name = clade_list[0].name
        sp = fstlib.align(lower_pass_fst, fsa_dict[normal_name], fsa_dict[root_name])
        fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

        logger.debug(f"Ancestor reconstruction using {n_cores} threads: Down the tree")
        # Down the tree (root to leaf), level by level
        for depth in sorted(levels.keys()):
            pairs = []
            for node in levels[depth]:
                if len(node.clades) != 0:
                    for child in node.clades:
                        if len(child.clades) != 0:
                            pairs.append((node.name, child.name))

            if not pairs:
                continue

            if len(pairs) > 1:
                results = Parallel(n_jobs=n_cores)(
                    delayed(_align_task)(fsa_dict[pname], fsa_dict[cname], lower_pass_fst)
                    for pname, cname in pairs
                )
                for (_,cname), result in zip(pairs, results):
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
    if not upper_cache:
        return fsa_dict
    else:
        return fsa_dict, uppass_cache

def reconstruct_ancestors_incremental(tree, samples_dict, upper_pass_fst, lower_pass_fst, normal_name, prune_weight, old_uppass_cache, nni_move):
    """
    incremental reconstruction of ancestors after a single NNI move.
    Only the affected dirty nodes' candidates (up the tree) fsas are reconstructed.
    """
    def __build_parent_of_and_clade_loop_up_table(tree, normal_name):
        parent_of = {}
        clade_look_up_table = {}
        for clade in tree.find_clades(order="preorder"):
            clade_look_up_table[clade.name] = clade
            for child in clade.clades:
                if child.name != normal_name:
                    parent_of[child.name] = clade.name
        return parent_of, clade_look_up_table

    def __dirty_node_list(parent_of, v_name, normal_name, clade_look_up_table):
        dirty = []
        node_name = v_name
        while node_name is not None and node_name != normal_name:
            dirty.append(clade_look_up_table[node_name])
            node_name = parent_of[node_name]
        return dirty

    # Get the dirty nodes (from v_name to root)
    parent_of, clade_look_up_table = __build_parent_of_and_clade_loop_up_table(tree, normal_name)
    dirty_nodes_l = __dirty_node_list(parent_of, nni_move.v_name, normal_name, clade_look_up_table)

    fsa_dict = old_uppass_cache | samples_dict
    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]

    logger.debug("Incremental ancestral reconstuction: Up the tree for only the dirty nodes")
    for node in dirty_nodes_l:
        children = [item for item in node.clades if item.name != normal_name]
        left_name = children[0].name
        right_name = children[1].name
        logger.debug(f"Clade: {node.name}, left: {left_name}, right: {right_name}")

        ## project
        intersection = intersect_clades_detmin(fsa_dict[left_name], fsa_dict[right_name], upper_pass_fst,
                                               prune_weight=prune_weight, detmin_before_intersect=False,
                                               detmin_after_intersect=True)
        fsa_dict[node.name] = intersection

    new_uppass_cache = {node.name: fsa_dict[node.name] for node in clade_list if len(node.clades) != 0}

    logger.debug("Incremental ancestral reconstuction for root")
    # root node is calculated separately w.r.t. normal node
    root_name = clade_list[0].name
    sp = fstlib.align(lower_pass_fst, fsa_dict[normal_name], fsa_dict[root_name])
    fsa_dict[root_name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    logger.info("Incremental ancestral reconstuction: Down the tree")
    # down the tree (root to leaf)
    for node in clade_list:
        if len(node.clades) != 0:
            children = [q for q in node.clades if len(q.clades) != 0]
            logger.debug(f"Clade: {node.name}, internal children: {children}")
            for child in children:
                sp = fstlib.align(lower_pass_fst, fsa_dict[node.name], fsa_dict[child.name])
                fsa_dict[child.name] = fstlib.arcmap(sp.copy().project('output'), map_type='rmweight')

    # check if ancestors were correctly reconstructed
    sample_lengths = {sample: len(medicc.tools.fsa_to_string(fsa_dict[sample])) for sample, fsa in fsa_dict.items()}
    normal_length = sample_lengths[normal_name]

    if np.any([x != normal_length for x in sample_lengths.values()]):
        raise MEDICCAncestorReconstructionError("Some ancestors could not be reconstructed. These are:\n"
                                                "{}".format('\n'.join([sample for sample, length in sample_lengths.items() if length != normal_length])) + \
                                                "\nCheck whether your normal sample contains segments with copy number zero")
    return fsa_dict, new_uppass_cache



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


