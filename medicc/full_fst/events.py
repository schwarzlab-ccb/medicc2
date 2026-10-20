"""One minimum-cost, replayable event history under the capped full-FST model.

The native scorer chooses the global WGD count and the intermediate profiles.
OpenFST ShortestPath then unfolds each small A or AD stage into L_LOH and G.
Loss/gain intervals may continue through zero-copy bins; gains never resurrect
an absent allele. WGD acts on the entire genome as min(2 * CN, 8).
"""
from dataclasses import dataclass
from functools import lru_cache
import ctypes as C
import json
import math

import pandas as pd
from fstlib.cext import pywrapfst as fst
from ._bindings import build_directory
from .pair_scorer import Model


@dataclass(frozen=True)
class Event:
    type: str
    start: int                 # Encoded genome column, inclusive; -1 for WGD.
    end: int                   # Encoded genome column, exclusive; -1 for WGD.
    stage: int                 # Number of preceding WGDs.


@dataclass
class History:
    cost: int
    wgds: int
    events: list[Event]


def replay(parent, events):
    """Apply ordered unit events, including zero-preserving gains and capped WGDs."""
    values = [None if c == 'X' else int(c) for c in parent]
    for event in events:
        if event.type == 'wgd':
            values = [None if v is None else min(2*v, 8) for v in values]
            continue
        if event.type not in ('loss', 'gain') or not 0 <= event.start < event.end <= len(values):
            raise ValueError('Invalid interval event')
        for i in range(event.start, event.end):
            if values[i] is None:
                raise ValueError('Event crosses a chromosome or allele boundary')
            if event.type == 'loss':
                values[i] = max(0, values[i]-1)
            elif values[i] > 0:
                values[i] += 1
                if values[i] > 8:
                    raise ValueError('Gain leaves the CN alphabet')
    return ''.join('X' if v is None else str(v) for v in values)


@lru_cache(maxsize=1)
def _parts():
    """Read the exact exported L_LOH/G primitives, not legacy MEDICC FST files."""
    lines = iter((build_directory('direct_pair')/'ancestral_parts.txt').read_text().splitlines())
    parts = {}
    for header in lines:
        name, size, start, count = header.split()
        size, start, count = int(size), int(start), int(count)
        finals = list(map(int, next(lines).split()))
        arcs = [[[] for _ in range(10)] for _ in range(size)]
        for _ in range(count):
            source, a, b, weight, target = map(int, next(lines).split())
            if name not in ('L_LOH', 'G'):
                continue
            # These primitive states count active unit intervals. Verify this
            # interpretation against EVERY arc instead of assuming state IDs.
            if a == 9:
                valid = b == 9 and target == 0 and weight == 0
            elif name == 'L_LOH':
                valid = b == max(0, a-target) and weight == max(0, target-source)
            else:
                valid = b == (a+target if a else 0) and weight == max(0, target-source)
            if not valid:
                raise ValueError('Unsupported event interpretation of ' + name)
            arcs[source][a].append((b, weight, target))
        if name in ('L_LOH', 'G'):
            if start != 0 or finals != [0]*size:
                raise ValueError('Unsupported primitive start/final weights')
            parts[name] = arcs
    return parts['L_LOH'], parts['G']


def _stage(parent, child, doubled):
    """OpenFST shortest path on a small, layered, endpoint-constrained L/G graph."""
    loss, gain = _parts()
    graph = fst.VectorFst()
    graph.add_state(); graph.set_start(0)
    current = {(0, 0): 0}
    weights = [fst.Weight('tropical', str(i)) for i in range(16)]
    for a, b in zip(parent, child):
        a, b = (9 if a == 'X' else int(a)), (9 if b == 'X' else int(b))
        following = {}
        for (ql, qg), source in current.items():
            for middle, wl, tl in loss[ql][a]:
                for output, wg, tg in gain[qg][middle]:
                    observed = min(2*output, 8) if doubled and output != 9 else output
                    if observed != b:
                        continue
                    key = (tl, tg)
                    if key not in following:
                        following[key] = graph.add_state()
                    # Nonzero labels encode the two active interval counts.
                    label = 1 + tl*8 + tg
                    graph.add_arc(source, fst.Arc(label, label, weights[wl+wg], following[key]))
        current = following
    for state in current.values():
        graph.set_final(state, weights[0])
    path = fst.shortestpath(graph)
    state = path.start()
    if state < 0:
        raise ValueError('No literal L_LOH/G decomposition of the chosen stage')
    labels, cost = [], 0
    while True:
        edges = list(path.arcs(state))
        if not edges:
            final = float(path.final(state))
            if not math.isfinite(final):
                raise ValueError('Event path is not accepting')
            cost += int(final)
            break
        if len(edges) != 1:
            raise ValueError('Expected one shortest event path')
        arc = edges[0]
        cost += int(float(arc.weight))
        if arc.ilabel:
            labels.append(divmod(arc.ilabel-1, 8))
        state = arc.nextstate
    if len(labels) != len(parent):
        raise ValueError('Event path changed the genome length')
    return labels, cost


def _intervals(counts, word, kind, stage):
    """Turn active interval counts into ordered unit events; separators close all."""
    opened, events = [], []
    for i in range(len(word)+1):
        count = counts[i] if i < len(word) and word[i] != 'X' else 0
        while len(opened) > count:
            events.append(Event(kind, opened.pop(), i, stage))
        while len(opened) < count:
            opened.append(i)
    return events


class EventCaller:
    """Own one exact pair scorer; event calling is optional and has no search cap."""
    def __init__(self):
        self.scorer = Model(backend='direct_pair')
        self.native = self.scorer.lib.am_pair_history
        self.native.argtypes = [C.c_void_p, C.c_char_p, C.c_char_p, C.c_int]
        self.native.restype = C.c_char_p

    def close(self):
        self.scorer.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def branch(self, parent, child, *, expected_score=None):
        if not parent or len(parent) != len(child):
            raise ValueError('Event calling requires equal, nonempty genomes')
        if any(c not in '012345678X' for c in parent+child):
            raise ValueError('Invalid CN symbol')
        if any((a == 'X') != (b == 'X') for a, b in zip(parent, child)):
            raise ValueError('Chromosome/allele boundary mismatch')
        # Identical adjacent endpoint pairs can share a path. Expand its states
        # back to ALL original bins before emitting genomic event coordinates.
        keep = [i for i in range(len(parent)) if i == 0 or parent[i] == 'X'
                or (parent[i], child[i]) != (parent[i-1], child[i-1])]
        x, y = ''.join(parent[i] for i in keep), ''.join(child[i] for i in keep)
        score = self.scorer.score(x, y, directed=True)
        if not math.isfinite(score):
            raise ValueError('No feasible directed capped event history')
        if expected_score is not None and score != expected_score:
            raise ValueError(f'Branch score {expected_score} differs from event model score {score}')
        if score == 0:
            if parent != child:
                raise ValueError('Nonidentity zero-cost branch')
            return History(0, 0, [])
        raw = self.native(self.scorer.h, x.encode('ascii'), y.encode('ascii'), int(score))
        if raw is None:
            raise RuntimeError(self.scorer.lib.am_error().decode())
        trace = json.loads(raw)
        if trace['cost'] != score:
            raise ValueError('Native event traceback differs from directed score')
        events, value, current = [], trace['wgds'], x
        for stage, target in enumerate(trace['stages']):
            counts, cost = _stage(current, target, stage < trace['wgds'])
            expanded = [None]*len(parent)
            for j, start in enumerate(keep):
                end = keep[j+1] if j+1 < len(keep) else len(parent)
                expanded[start:end] = [counts[j]]*(end-start)
            events.extend(_intervals([c[0] for c in expanded], parent, 'loss', stage))
            events.extend(_intervals([c[1] for c in expanded], parent, 'gain', stage))
            if stage < trace['wgds']:
                events.append(Event('wgd', -1, -1, stage))
            current = target
            value += cost
        if value != score or len(events) != score or replay(parent, events) != child:
            raise ValueError('Event replay or event count disagrees with the capped score')
        return History(int(score), trace['wgds'], events)


def call_tree_events(tree, profiles, coordinates, alleles):
    """Return a legacy-shaped event table with explicit parent, order and stage."""
    positions = []
    for allele in alleles:
        previous = None
        for coordinate in coordinates:
            chrom, start, end = coordinate
            if previous is not None and chrom != previous:
                positions.append(None)
            positions.append((allele, chrom, int(start), int(end)))
            previous = chrom
        positions.append(None)
    positions.pop()  # No separator after the final allele.
    records = []
    with EventCaller() as caller:
        for parent in tree.find_clades():
            for child in parent.clades:
                x, y = profiles[parent.name], profiles[child.name]
                if len(x) != len(positions) or any((c == 'X') != (p is None) for c, p in zip(x, positions)):
                    raise ValueError('Event coordinates do not match the encoded genome')
                history = caller.branch(x, y, expected_score=child.branch_length)
                for order, event in enumerate(history.events):
                    if event.type == 'wgd':
                        allele, chrom = 'both', 'chr0'
                        start = min(c[1] for c in coordinates)
                        end = max(c[2] for c in coordinates)
                    else:
                        allele, chrom, start, _ = positions[event.start]
                        last_allele, last_chrom, _, end = positions[event.end-1]
                        if (allele, chrom) != (last_allele, last_chrom):
                            raise ValueError('Event crosses an allele/chromosome boundary')
                    records.append(dict(sample_id=child.name, chrom=chrom, start=start, end=end,
                                        allele=allele, type=event.type,
                                        cn_child=0 if event.type == 'wgd' else int(y[event.start]),
                                        parent_id=parent.name,
                                        event_order=order+1, stage=event.stage, cost=1,
                                        genome_start=event.start, genome_end=event.end))
    columns = ['sample_id', 'chrom', 'start', 'end', 'allele', 'type', 'cn_child', 'parent_id',
               'event_order', 'stage', 'cost', 'genome_start', 'genome_end']
    return pd.DataFrame(records, columns=columns).set_index(['sample_id', 'chrom', 'start', 'end'])
