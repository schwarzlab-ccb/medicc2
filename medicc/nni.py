"""Nearest-neighbour interchange (NNI) refinement of a tree topology.

An NNI move at an internal edge (u, v), where v is an internal child of u, swaps the
other child of u with one of the two children of v. Every such edge below the root
gives two moves. The normal sample stays the outgroup and no internal node is created,
renamed or removed.

The search is independent of the copy-number model: it takes a function that scores
a tree (the total tree length after ancestor reconstruction) and runs a steepest
descent over NNI neighbours.
"""
import copy
import logging

logger = logging.getLogger(__name__)


def nni_neighbours(tree, normal_name='diploid'):
    """Generate all NNI neighbours of a rooted binary tree.

    For every internal node u with exactly two children other than the normal, and every
    internal child v of u, the other child of u is swapped with each child of v in turn.
    Moves are generated in preorder of u, then in child order. The input tree is not changed.

    Parameters
    ----------
    tree : Bio.Phylo tree
        Rooted tree, e.g. the tree used for ancestor reconstruction, whose root holds the
        normal and the two subtrees.
    normal_name : str
        Name of the normal sample, which is never moved.

    Yields
    ------
    Bio.Phylo tree
        A modified deep copy of ``tree`` for each move.
    """
    for u_index, u in enumerate(tree.find_clades()):
        children = [i for i, child in enumerate(u.clades) if child.name != normal_name]
        if len(children) != 2:
            continue
        for v_index in children:
            v = u.clades[v_index]
            if len(v.clades) < 2:
                continue
            other_index = children[1] if v_index == children[0] else children[0]
            for target_index in range(len(v.clades)):
                neighbour = copy.deepcopy(tree)
                new_u = list(neighbour.find_clades())[u_index]
                new_v = new_u.clades[v_index]
                new_u.clades[other_index], new_v.clades[target_index] = (
                    new_v.clades[target_index], new_u.clades[other_index])
                yield neighbour


def topology_key(tree):
    """Hashable key of the rooted topology, independent of internal node names and child order."""
    def key(clade):
        if clade.is_terminal():
            return (0, clade.name)
        return (1, tuple(sorted(key(child) for child in clade.clades)))

    return key(tree.root)


def nni_search(tree, evaluate, max_iter=None, normal_name='diploid'):
    """Steepest-descent NNI search that minimises the score of a tree.

    Each sweep evaluates every NNI neighbour (see ``nni_neighbours``) of the current tree
    whose topology was not evaluated before and moves to the best neighbour if it is
    strictly better than the current tree. Ties are broken by enumeration order, so the
    result is deterministic. The search stops at a local optimum or when ``max_iter``
    candidates have been evaluated; the best candidate found until then is kept.

    Parameters
    ----------
    tree : Bio.Phylo tree
        Start tree (not changed).
    evaluate : callable
        Function tree -> score (float), e.g. the total tree length after ancestor
        reconstruction. It must not change the tree it gets.
    max_iter : int, optional
        Maximum number of evaluated candidate topologies (the start tree is not counted).
        Unlimited by default.
    normal_name : str
        Name of the normal sample.

    Returns
    -------
    best_tree : Bio.Phylo tree
        The tree with the lowest score found.
    trace : list of dict
        One entry per accepted step with the keys ``step`` (0 is the start tree),
        ``tree_length`` (score) and ``evaluated`` (candidates evaluated so far).

    Raises
    ------
    ValueError
        If ``max_iter`` is smaller than 1.
    """
    if max_iter is not None and max_iter < 1:
        raise ValueError("max_iter must be at least 1.")

    best_tree = tree
    best_score = evaluate(tree)
    seen = {topology_key(tree)}
    evaluated = 0
    limit_reached = False
    trace = [dict(step=0, tree_length=best_score, evaluated=0)]
    logger.info(f"NNI start: tree length {best_score}.")

    while not limit_reached:
        winner, winner_score = None, best_score
        for neighbour in nni_neighbours(best_tree, normal_name):
            key = topology_key(neighbour)
            if key in seen:
                continue
            if max_iter is not None and evaluated >= max_iter:
                limit_reached = True
                break
            seen.add(key)
            score = evaluate(neighbour)
            evaluated += 1
            logger.debug(f"NNI candidate {evaluated}: tree length {score}")
            if score < winner_score:
                winner, winner_score = neighbour, score
        if winner is None:
            break
        best_tree, best_score = winner, winner_score
        trace.append(dict(step=len(trace), tree_length=best_score, evaluated=evaluated))
        logger.info(f"NNI step {len(trace) - 1}: tree length {best_score} "
                    f"({evaluated} candidates evaluated).")

    if limit_reached:
        logger.info(f"NNI stopped at the limit of {max_iter} evaluated candidates.")
    else:
        logger.info(f"NNI reached a local optimum after {evaluated} evaluated candidates.")
    return best_tree, trace
