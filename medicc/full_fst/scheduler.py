"""Bottom-up message joins, then top-down profile choices on a fixed binary tree.

Each worker has its own native Model. Completed child messages are shared read-only.
Ready nodes start as soon as their dependencies finish. The original level
scheduler remains available for controlled comparisons.
"""

import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed, wait, FIRST_COMPLETED
from .model import Model, ExecutionBudget


def reconstruct(
    children,
    leaves,
    root,
    normal,
    workers=64,
    state_limit=0,
    retain_messages=False,
    progress=None,
    model_factory=None,
    scheduling="ready",
    branch_parallel=False,
):
    if scheduling not in {"ready", "levels"}:
        raise ValueError("Expected ready or levels scheduling")
    if workers < 1:
        raise ValueError("workers must be positive")
    template = lambda s: "".join("X" if c == "X" else "#" for c in s)
    if len({template(s) for s in [normal, *leaves.values()]}) != 1:
        raise ValueError("Column/separator mismatch")
    levels = {}
    seen = set()

    def walk(v, d):
        if v in seen:
            raise ValueError("Expected a tree")
        seen.add(v)
        if v in children:
            if v in leaves or len(children[v]) != 2:
                raise ValueError("Expected binary tree")
            levels.setdefault(d, []).append(v)
            for c in children[v]:
                walk(c, d + 1)
        elif v not in leaves:
            raise ValueError("Missing leaf")

    walk(root, 0)
    if seen != set(children) | set(leaves):
        raise ValueError("Unreachable nodes")
    budget = ExecutionBudget(workers) if branch_parallel else None
    local = threading.local()
    models = []
    messages = {}
    profiles = dict(leaves)
    counts = dict(up=0, down=0)

    def model():
        if not hasattr(local, "model"):
            local.model = (model_factory or Model)(state_limit=state_limit)
            local.model.execution_budget = budget
            models.append(local.model)
        return local.model

    def emit(**kw):
        if progress:
            progress(kw)

    def up(v):
        m = model()
        start = time.perf_counter()
        a, b = children[v]
        x = m.string(leaves[a]) if a in leaves else messages[a]
        y = m.string(leaves[b]) if b in leaves else messages[b]
        try:
            out = m.join(x, y)
        except Exception as e:
            raise RuntimeError(f"up {v}: {e}; stats={m.stats()}") from e
        return v, out, time.perf_counter() - start, m.stats()

    def down(v, parent):
        m = model()
        start = time.perf_counter()
        try:
            s, cost = m.resolve(parent, messages[v])
        except Exception as e:
            raise RuntimeError(f"down {v}: {e}; stats={m.stats()}") from e
        return v, s, time.perf_counter() - start, cost, m.stats()

    parents = {c: v for v, cs in children.items() for c in cs}
    depths = {v: depth for depth, nodes in levels.items() for v in nodes}
    durations = {"up": {}, "down": {}}
    peak_retained = 0

    def accept_up(future):
        nonlocal peak_retained
        v, msg, sec, stats = future.result()
        messages[v] = msg
        peak_retained += msg.size
        counts["up"] += 1
        durations["up"][v] = sec
        emit(
            phase="up",
            node=v,
            seconds=sec,
            states=msg.size,
            completed=counts["up"],
            depth=depths[v],
            stats=stats,
        )
        return v

    def accept_down(future):
        v, s, sec, cost, stats = future.result()
        profiles[v] = s
        counts["down"] += 1
        durations["down"][v] = sec
        emit(
            phase="down",
            node=v,
            seconds=sec,
            score=cost,
            completed=counts["down"],
            depth=depths[v],
            stats=stats,
        )
        if not retain_messages:
            messages.pop(v).close()
        return v

    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="ancestor") as pool:
        start = time.perf_counter()
        if scheduling == "levels":
            for depth in sorted(levels, reverse=True):
                jobs = [pool.submit(up, v) for v in levels[depth]]
                for future in as_completed(jobs):
                    accept_up(future)
        else:
            # The coordinator alone publishes results and updates dependency counts.
            # Tasks only read child messages that have already been published.
            pending_children = {
                v: sum(c in children for c in cs) for v, cs in children.items()
            }
            pending = {
                pool.submit(up, v)
                for depth in sorted(levels, reverse=True)
                for v in levels[depth]
                if pending_children[v] == 0
            }
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    v = accept_up(future)
                    parent = parents.get(v)
                    if parent is not None:
                        pending_children[parent] -= 1
                        if pending_children[parent] == 0:
                            pending.add(pool.submit(up, parent))
        up_seconds = time.perf_counter() - start
        start = time.perf_counter()
        if scheduling == "levels":
            for depth in sorted(levels):
                jobs = [
                    pool.submit(down, v, normal if v == root else profiles[parents[v]])
                    for v in levels[depth]
                ]
                for future in as_completed(jobs):
                    accept_down(future)
        else:
            pending = {pool.submit(down, root, normal)} if root in children else set()
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    v = accept_down(future)
                    for child in children[v]:
                        if child in children:
                            pending.add(pool.submit(down, child, profiles[v]))
        down_seconds = time.perf_counter() - start

    # These are observed task durations, not CPU utilization measurements. Their
    # dependency-chain sums give a lower bound on each pass's ideal schedule.
    critical = {}
    for phase in ("up", "down"):
        chains = {}
        for depth in sorted(levels, reverse=True):
            for v in levels[depth]:
                chains[v] = durations[phase][v] + max(
                    (chains.get(c, 0) for c in children[v]), default=0
                )
        critical[phase] = chains.get(root, 0)
    metrics = dict(
        up_seconds=up_seconds,
        down_seconds=down_seconds,
        retained_internal_nodes_peak=peak_retained,
        workers_requested=workers,
        worker_models=len(models),
        up_nodes=counts["up"],
        down_nodes=counts["down"],
        scheduling=scheduling,
        task_seconds={
            phase: sum(values.values()) for phase, values in durations.items()
        },
        critical_path_seconds=critical,
        branch_parallel=branch_parallel,
        peak_native_workers=budget.peak if budget else None,
    )
    # Messages keep their originating model alive when requested by validation.
    if not retain_messages:
        for m in models:
            m.close()
    return profiles, messages, metrics
