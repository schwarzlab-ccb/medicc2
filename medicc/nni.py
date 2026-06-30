"""NNI (Nearest Neighbor Interchange) move enumeration for MEDICC2.

Mirrors the shape and conventions of medicc/spr.py. Yields all NNI neighbors of
a rooted MEDICC tree along with a synthetic SPRResult for each, so they can be
fed into reconstruct_ancestors_incremental without any new infrastructure.

The NNI move at internal edge (u, v) — where u is the parent of v and v is itself
internal — picks the "uncle" A (the other child of u) and swaps it with one of v's
children (B or C). Two distinct moves per qualifying edge.

Edges incident to the root are excluded, matching the SPR convention in spr.py.
The MEDICC search-shape tree has the diploid as root with a single structural
child (the MRCA), so root and root-children are skipped.
"""
import copy

from medicc.spr import SPRResult

# ---------------------------------------------------------------------------
# Move descriptor — a lightweight picklable description of one NNI move.
# Used to defer deepcopy into worker processes for parallel evaluation.
# ---------------------------------------------------------------------------

class _NNIMove:
    """Identifies one NNI move by the names of the four involved nodes."""
    __slots__ = ("u_name", "v_name", "A_name", "swap_target_name")

    def __init__(self, u_name, v_name, A_name, swap_target_name):
        self.u_name = u_name
        self.v_name = v_name
        self.A_name = A_name
        self.swap_target_name = swap_target_name


def _is_root_or_root_child(clade, root, root_children):
    return clade is root or clade in root_children


def _apply_nni_swap(tree, u_name, v_name, A_name, swap_target_name):
    """Deepcopy `tree` and apply an NNI swap at edge (u, v), exchanging A and swap_target.

    Pre: u is parent of v; A is u's other child; swap_target is one of v's children.
    Post: u's children are {swap_target, v}; v's children are {A, the_other_v_child}.

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
    """Return a list of _NNIMove descriptors for all valid NNI moves on `tree`.

    Does NOT deepcopy the tree — just collects node names.  Used by both the
    serial path and the parallel path to separate enumeration from execution.
    """
    root = tree.root
    root_children = set(root.clades)
    moves = []
    for u in tree.find_clades():
        if _is_root_or_root_child(u, root, root_children):
            continue
        if len(u.clades) < 2:
            continue
        for v in u.clades:
            if len(v.clades) < 2:
                continue
            siblings_of_v = [c for c in u.clades if c is not v]
            if not siblings_of_v:
                continue
            A = siblings_of_v[0]
            for swap_target in list(v.clades):
                moves.append(_NNIMove(u.name, v.name, A.name, swap_target.name))
    return moves


def nni_neighbors(tree):
    """Yield (new_tree, synthetic_SPRResult) for every NNI neighbor of `tree`.

    Iterates over internal edges (u, v) skipping edges incident to the root,
    and for each qualifying edge yields the 2 uncle-swap variants.

    The synthetic SPRResult uses prune_grandparent=u.name and regraft_parent=v.name
    so reconstruct_ancestors_incremental._get_dirty_nodes_up walks
    {v, u, parent(u), ..., root} — exactly the dirty set for an NNI move.
    """
    for move in _enumerate_moves(tree):
        new_tree = _apply_nni_swap(
            tree=tree,
            u_name=move.u_name,
            v_name=move.v_name,
            A_name=move.A_name,
            swap_target_name=move.swap_target_name,
        )
        synthetic = SPRResult(
            tree=new_tree,
            pruned_subtree_root=move.swap_target_name,  # cosmetic only
            prune_grandparent=move.u_name,
            regraft_parent=move.v_name,
        )
        yield new_tree, synthetic


def _eval_nni_neighbor(tree, move, old_uppass_cache, samples_dict, fst,
                       normal_name, prune_weight, step=0,
                       upper_pass_fst=None, lower_pass_fst=None):
    """Apply one NNI move and evaluate it.  Runs inside a worker process.

    Deepcopy happens here so the main process never copies the tree before
    forking.  Returns (new_tree, ancestors, uppass_cache, score, step).

    upper_pass_fst/lower_pass_fst: optional separate FSTs for the up-pass and
    down-pass (the "lower-half" length-encoding objective). When None both fall
    back to `fst`. Branch lengths / the move score are computed under the
    down-pass FST so the hill-climb ranks neighbors by the length-encoding
    objective (e.g. open=5000: 5000*events + span).
    """
    # Import lazily — these are heavy and only needed inside workers.
    from medicc.ancestors import reconstruct_ancestors_incremental
    import medicc.tools
    import medicc.core as core

    score_fst = lower_pass_fst if lower_pass_fst is not None else fst

    new_tree = _apply_nni_swap(
        tree=tree,
        u_name=move.u_name,
        v_name=move.v_name,
        A_name=move.A_name,
        swap_target_name=move.swap_target_name,
    )
    synthetic = SPRResult(
        tree=new_tree,
        pruned_subtree_root=move.swap_target_name,
        prune_grandparent=move.u_name,
        regraft_parent=move.v_name,
    )
    ancestors, new_uppass_cache = reconstruct_ancestors_incremental(
        tree=new_tree,
        old_uppass_cache=old_uppass_cache,
        samples_dict=samples_dict,
        fst=fst,
        normal_name=normal_name,
        spr_result=synthetic,
        prune_weight=prune_weight,
        upper_pass_fst=upper_pass_fst,
        lower_pass_fst=lower_pass_fst,
    )
    core.update_branch_lengths(new_tree, score_fst, ancestors, normal_name)
    score = medicc.tools.sum_of_branch_length(new_tree)
    return new_tree, ancestors, new_uppass_cache, score, step


def evaluate_nni_neighbors_parallel(tree, old_uppass_cache, samples_dict, fst,
                                    normal_name, prune_weight, n_cores,
                                    step_start=0, upper_pass_fst=None,
                                    lower_pass_fst=None):
    """Evaluate all NNI neighbors in parallel using joblib.

    Returns a list of (new_tree, ancestors, uppass_cache, score, step) tuples,
    one per NNI neighbor.  The deepcopy of `tree` is deferred into each
    worker so the main process does not copy the tree before forking.
    `step_start` sets the step number of the first neighbor; subsequent
    neighbors get step_start+1, step_start+2, etc.

    upper_pass_fst/lower_pass_fst: optional separate up-/down-pass FSTs for the
    length-encoding objective (see _eval_nni_neighbor). Forwarded to each worker.
    """
    try:
        from joblib import Parallel, delayed
    except ImportError:
        raise ImportError("joblib must be installed for parallel NNI evaluation")

    moves = _enumerate_moves(tree)
    if not moves:
        return []

    results = Parallel(n_jobs=n_cores)(
        delayed(_eval_nni_neighbor)(
            tree, move, old_uppass_cache, samples_dict, fst, normal_name, prune_weight,
            step=step_start + i,
            upper_pass_fst=upper_pass_fst,
            lower_pass_fst=lower_pass_fst,
        )
        for i, move in enumerate(moves)
    )
    return results
