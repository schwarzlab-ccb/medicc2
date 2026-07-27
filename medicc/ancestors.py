import logging

import Bio
import Bio.Phylo
import fstlib
import numpy as np

import medicc

logger = logging.getLogger(__name__)

def reconstruct_ancestors(tree, samples_dict, upper_pass_fst, lower_pass_fst, normal_name, prune_weight=0):
    '''
    Reconstruct ancestors using upper fst for traversing up the tree and lower fst to traverse down the tree.
    '''

    if len(samples_dict) == 2:
        return samples_dict

    fsa_dict = samples_dict.copy()
    tree = Bio.Phylo.BaseTree.copy.deepcopy(tree)

    clade_list = [clade for clade in tree.find_clades(order="preorder") if clade.name != normal_name]
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

    # check if ancestors were correctly reconstructed
    sample_lengths = {sample: len(medicc.tools.fsa_to_string(fsa_dict[sample])) for sample, fsa in fsa_dict.items()}
    normal_length = sample_lengths[normal_name]

    if np.any([x != normal_length for x in sample_lengths.values()]):
        raise MEDICCAncestorReconstructionError("Some ancestors could not be reconstructed. These are:\n"
                                                "{}".format('\n'.join([sample for sample, length in sample_lengths.items() if length != normal_length])) + \
                                                "\nCheck whether your normal sample contains segments with copy number zero")

    return fsa_dict


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
