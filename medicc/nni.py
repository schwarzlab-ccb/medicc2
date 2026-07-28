"""
NNI (nearest neighbor interchange) move for MEDICC2

The NNI move at internal edge (u, v), where u is the parent of v and v is itself
internal - picksthe "uncle" A (the other child of u) and swaps it with one of v's
children (B or C). Two distinct moves per qualifying edge.
"""

import copy

class _NNIMove:
    """Identifies one NNI move by the names of the four involved nodes"""
    __slots__ = ("u_name", "v_name", "A_name", "swap_target_name")

    def __init__(self, u_name, v_name, A_name, swap_target_name):
        self.u_name = u_name
        self.v_name = v_name
        self.A_name = A_name
        self.swap_target_name = swap_target_name

def _is_root(clade, root):
    return clade is root

def _apply_nni_swap(tree, u_name, v_name, A_name, swap_target_name):
    """
    Deepcopy 'tree' and apply an NNI swap at edge (u, v), exchanging A and swap_targe.

    Pre: u is parent of v; A is u's other child; swap_target is one of v's children.
    post: u's children are {swap_target, v}; v's children are {A, the_other_v_child}.

    Returns the new tree. No new internal nodes are created; all names preserved.
    """

    new_tree = copy.deepcopy(tree)

    # Locate u, v, A, swap_target by name in the copy
    name_to_clade = {c.name: c for c in new_tree.find_clades() if c.name is not None}
    u = name_to_clade[u_name]
    v = name_to_clade[v_name]
    A = name_to_clade[A_name]
    swap_target = name_to_clade[swap_target_name]

    # Rewire: remove A from u, append swap_target to u
    u.clades = [c for c in u.clades if c is not A]
    u.clades.append(swap_target)

    # Rewire: remove swap_target from v, append A to v
    v.clades = [c for c in v.clades if c is not swap_target]
    v.clades.append(A)

    return new_tree

def _enumerate_moves(tree):
    """
    Return a list of _NNIMove descriptors for all valid NNI moves on `tree`.
    """

    root = tree.root
    root_children = set(root.clades)
    moves = []

    # Depth first traversal
    for u in tree.find_clades():
        if _is_root(u, root):
            continue
        if len(u.clades) < 2:
            continue
        for v in u.clades:
            if len(v.clades) < 2:
                continue
            siblings_of_v = [c for c in u.clades if c is not v]
            A = siblings_of_v[0]
            for swap_target in list(v.clades):
                moves.append(_NNIMove(u.name, v.name, A.name, swap_target.name))
    return moves


def nni_neighbors(tree):
    """Yield (new_tree, _NNIMove) for each valid NNI move on 'tree'."""

    for move in _enumerate_moves(tree):
        new_tree = _apply_nni_swap(
            tree=tree,
            u_name=move.u_name,
            v_name=move.v_name,
            A_name=move.A_name,
            swap_target_name=move.swap_target_name
        )
        yield new_tree, move

def _eval_nni_neighbor(tree, move, old_uppass_cache, samples_dict, upper_pass_fst, lower_pass_fst,
                       normal_name, prune_weight, visited_solutions, step=0):
    """
    Apply one NNI move and evaluate it. Runs inside a worker process.

    Returns (new_tree, ancestors, uppass_cache, score, step, nni_move)
    """
    # Lazy import
    import medicc.ancestors
    import medicc.tools
    import medicc.core

    new_tree = _apply_nni_swap(tree, move.u_name, move.v_name, move.A_name, move.swap_target_name)
    new_tree_hash = medicc.tree_hash.get_topology_hash(new_tree)
    if new_tree_hash in visited_solutions:
        score = visited_solutions[new_tree_hash]
        return new_tree, None, None, score, step, move
    else:
        ancestors, new_uppass_cache = medicc.ancestors.reconstruct_ancestors_incremental(
            tree=new_tree,
            samples_dict=samples_dict,
            upper_pass_fst=upper_pass_fst,
            lower_pass_fst=lower_pass_fst,
            normal_name=normal_name,
            prune_weight=prune_weight,
            old_uppass_cache=old_uppass_cache,
            nni_move=move,
        )
        medicc.core.update_branch_lengths(new_tree, lower_pass_fst, ancestors, normal_name)
        score = medicc.tools.sum_of_branch_length(new_tree)
        return new_tree, ancestors, new_uppass_cache, score, step, move


def evaluate_nni_neighbors_parallel(tree, old_uppass_cache, samples_dict,
                                    normal_name, prune_weight, n_cores, upper_pass_fst, lower_pass_fst, visited_solutions,
                                    step_start=0):
    """
    Evaluate all NNI neighbors in parallel using joblib

    Returns a list of (new_tree, ancestors, uppass_cache, score, step) tuples,
    one per NNI neighbor.
    """

    from joblib import Parallel, delayed

    moves = _enumerate_moves(tree)

    results = Parallel(n_jobs=n_cores)(
        delayed(_eval_nni_neighbor)(
            tree, move, old_uppass_cache, samples_dict, upper_pass_fst, lower_pass_fst,
            normal_name, prune_weight, visited_solutions, step=step_start+i,
        )
        for i, move in enumerate(moves)
    )
    return results


