"""Threaded pair distances and ancestor reconstruction with the native model.

Every worker thread owns its own native Model. The coordinator only hands out
work whose inputs are complete; messages are immutable once published.
"""
import logging
import threading
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait

import numpy as np
import pandas as pd

from medicc.native import shortening
from medicc.native.model import Model

logger = logging.getLogger(__name__)


class _WorkerModels:
    """Create one Model per worker thread on first use; close all of them at the end."""

    def __init__(self, config):
        self.config = config
        self.local = threading.local()
        self.lock = threading.Lock()
        self.models = []

    def get(self):
        if not hasattr(self.local, 'model'):
            self.local.model = Model(self.config)
            with self.lock:
                self.models.append(self.local.model)
        return self.local.model

    def close(self):
        for model in self.models:
            model.close()


def pairwise_distances(cn_str_dict, config, n_cores=1):
    """Calculate the symmetric native distances between all pairs of profiles.

    Args:
        cn_str_dict (dict): Sample name -> copy-number profile (chromosomes separated by X).
        config (ModelConfig): Native model selection.
        n_cores (int, optional): Number of worker threads. Defaults to 1.

    Returns:
        pandas.DataFrame: Distance matrix indexed by the sample names in input order.
        Pairs that cannot be related have distance inf.
    """
    samples = list(cn_str_dict.keys())
    n = len(samples)
    distances = np.zeros((n, n), dtype=float)
    prepared = [shortening.prepare(cn_str_dict[sample]) for sample in samples]
    models = _WorkerModels(config)

    def score_row(i):
        model = models.get()
        row = [model.score(*shortening.pair(prepared[i], prepared[j])) for j in range(i + 1, n)]
        model.clear()  # caches are only reused within a row
        return i, row

    try:
        with ThreadPoolExecutor(max_workers=n_cores or 1) as pool:
            for future in [pool.submit(score_row, i) for i in range(n - 1)]:
                i, row = future.result()
                distances[i, i + 1:] = row
                distances[i + 1:, i] = row
    finally:
        models.close()

    return pd.DataFrame(distances, index=samples, columns=samples)


def reconstruct(children, leaves, root, normal, config, n_cores=1):
    """Reconstruct one optimal set of ancestral profiles on a rooted binary tree.

    The up pass joins the child messages of every internal node (all minimum-cost
    ancestors), the down pass resolves the root against the normal and then every
    internal node against its resolved parent.

    Args:
        children (dict): Internal node name -> list of its two child names (normal excluded).
        leaves (dict): Leaf name -> copy-number profile.
        root (str): Name of the root (the most recent common ancestor of all leaves).
        normal (str): Copy-number profile of the normal sample (parent of the root).
        config (ModelConfig): Native model selection.
        n_cores (int, optional): Number of worker threads. Defaults to 1.

    Returns:
        dict: Profiles of all leaves and internal nodes.

    Raises:
        ValueError: If the tree is not binary or its leaves do not match ``leaves``.
        RuntimeError: If the native engine fails, e.g. no common ancestor exists.
    """
    parents = {child: node for node, cs in children.items() for child in cs}
    for node, cs in children.items():
        if len(cs) != 2 or node in leaves:
            raise ValueError("Ancestor reconstruction requires a binary tree")
    if (set(parents) | {root}) != (set(children) | set(leaves)) or root in parents:
        raise ValueError("The tree does not match the samples")

    models = _WorkerModels(config)
    messages = {}
    profiles = dict(leaves)

    def join(node):
        model = models.get()
        temporary, operands = [], []
        for child in children[node]:
            if child in leaves:
                temporary.append(model.string(leaves[child]))
                operands.append(temporary[-1])
            else:
                operands.append(messages[child])
        try:
            return node, model.join(*operands)
        finally:
            for message in temporary:
                message.close()

    def resolve(node, parent_profile):
        return node, models.get().resolve(parent_profile, messages[node])[0]

    try:
        with ThreadPoolExecutor(max_workers=n_cores or 1, thread_name_prefix='ancestor') as pool:
            logger.info("Ancestor reconstruction: Up the tree")
            waiting = {node: sum(c in children for c in cs) for node, cs in children.items()}
            pending = {pool.submit(join, node) for node in children if waiting[node] == 0}
            while pending:
                finished, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in finished:
                    node, message = future.result()
                    messages[node] = message
                    parent = parents.get(node)
                    if parent is not None:
                        waiting[parent] -= 1
                        if waiting[parent] == 0:
                            pending.add(pool.submit(join, parent))

            logger.info("Ancestor reconstruction: Down the tree")
            pending = {pool.submit(resolve, root, normal)} if root in children else set()
            while pending:
                finished, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in finished:
                    node, profile = future.result()
                    profiles[node] = profile
                    messages.pop(node).close()
                    pending.update(pool.submit(resolve, child, profile)
                                   for child in children[node] if child in children)
    finally:
        for message in messages.values():
            message.close()
        models.close()

    return profiles
