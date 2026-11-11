"""Literal-FST test oracle, built with fstlib (OpenFST) from MEDICC's factory primitives.

One directed edge is the chain A (W1 A)^k with A = L_LOH * G and W1 = one optional
capped WGD (c -> min(2c, 8), cost 1). So T(x, y) is the cheapest history from x to
y with at most k WGDs. OpenFST is the reference; nothing here shares code with the
native engines or with scripts/generate_capped_ad.py.

The chain is never composed into one transducer (A W1 A W1 A is far too large).
Instead a weighted word language is pushed through one primitive at a time, and
determinized after each step, so only the words reachable from the input remain.
"""
from collections.abc import Mapping, Sequence
from functools import lru_cache
import math

import fstlib
import medicc

MAX_CN = 8
SEPARATOR = 'X'


def get_parts():
    """Return (primitives, symbol table): MEDICC's L_LOH, G, ... plus W1."""
    sym = medicc.create_symbol_table(MAX_CN, SEPARATOR)
    p = medicc.create_copynumber_fst(symbol_table=sym, sep=SEPARATOR, enable_wgd=True,
                                     max_num_wgds=3, wgd_x2=True, output_all=True)
    p['W1'] = medicc.factory.create_1step_WGD_fst(sym, SEPARATOR, wgd_cost=1, minimize=True, wgd_x2=True, wgd_cap=True)
    # Keep L_LOH and G separate. This independently checks the detmin NOWGD path.
    return p, sym


def doubling(sym, cost=0):
    """Mandatory capped WGD: every symbol c -> min(2c, MAX_CN), X -> X, `cost` once per genome."""
    f = fstlib.Fst()
    f.set_input_symbols(sym)
    f.set_output_symbols(sym)
    first, rest = f.add_state(), f.add_state()
    f.set_start(first)
    f.set_final(rest, 0)
    for c in [str(i) for i in range(MAX_CN + 1)] + [SEPARATOR]:
        d = c if c == SEPARATOR else str(min(2 * int(c), MAX_CN))
        f.add_arc(first, (c, d, cost, rest))
        f.add_arc(rest, (c, d, 0, rest))
    return f


def chain(parts, wgd, wgds):
    """The edge model A (wgd A)^wgds as a list of primitives (L_LOH, G, wgd, L_LOH, G, ...).

    `wgd` is parts['W1'] (optional WGD, cost 1) or doubling(sym) (mandatory WGD).
    """
    a = [parts['L_LOH'], parts['G']]
    return a + [wgd, *a] * wgds


class Oracle:
    """Exact scores and ancestor languages for short words, computed with OpenFST."""

    def __init__(self, wgds=2):
        self.parts, self.sym = get_parts()
        self.model = chain(self.parts, self.parts['W1'], wgds)

    def language(self, words: Mapping[str, float]) -> fstlib.Fst:
        """Weighted acceptor: each word with its cost."""
        out = None
        for word, cost in words.items():
            f = fstlib.factory.from_string(word, final_weight=cost, isymbols=self.sym, osymbols=self.sym)
            out = f if out is None else out.union(f)
        return fstlib.determinize(out.rmepsilon())

    def children(self, parents: fstlib.Fst, model: Sequence[fstlib.Fst] = None) -> fstlib.Fst:
        """Every child word y with cost min_x parents[x] + T(x, y)."""
        for p in model or self.model:
            parents = _step(parents, p, 'output')
        return parents

    def ancestors(self, children: fstlib.Fst, model: Sequence[fstlib.Fst] = None) -> fstlib.Fst:
        """Every ancestor word z with cost min_y T(z, y) + children[y]."""
        for p in reversed(model or self.model):
            children = _step(children, p, 'input')
        return children

    def score(self, x: str, y: str, model: Sequence[fstlib.Fst] = None) -> float:
        """Cheapest path from x to y through `model` (default: the edge model)."""
        return _total(fstlib.intersect(self.children(self.language({x: 0}), model).arcsort('olabel'),
                                       self.language({y: 0})))

    @lru_cache(maxsize=None)
    def directed(self, x: str, y: str) -> float:
        """T(x, y): cheapest edge history from parent x to child y."""
        return self.score(x, y)

    @lru_cache(maxsize=None)
    def kernel(self, x: str, y: str) -> float:
        """K(x, y) = min_z T(z, x) + T(z, y): the pair distance."""
        return self.pair(x, y)

    def pair(self, x: str, y: str, left: Sequence[fstlib.Fst] = None, right: Sequence[fstlib.Fst] = None) -> float:
        """min_z left(z, x) + right(z, y). Both default to the edge model."""
        zx = self.ancestors(self.language({x: 0}), left)
        zy = self.ancestors(self.language({y: 0}), right)
        return _total(fstlib.intersect(zx.arcsort('olabel'), zy))

    def optimum(self, a: Mapping[str, float], b: Mapping[str, float]) -> dict[str, int]:
        """All optimal ancestors z of two weighted child languages, with their costs.

        cost(z) = min_w (T(z, w) + a[w]) + min_w (T(z, w) + b[w]); only the minimum-cost z are kept.
        """
        joint = fstlib.intersect(self.ancestors(self.language(a)).arcsort('olabel'),
                                 self.ancestors(self.language(b)))
        if not math.isfinite(_total(joint)):
            return {}
        best = fstlib.determinize(joint)
        best.prune(weight=0)  # keep every word that ties the optimum
        return {word: int(cost) for word, cost in _words(best, self.sym).items()}


def _step(words, primitive, side):
    """Push a weighted word language through one primitive.

    side='output': words are inputs of `primitive`; keep its outputs (parent -> child).
    side='input':  words are outputs of `primitive`; keep its inputs (child -> parent).
    """
    if side == 'output':
        f = fstlib.compose(words.arcsort('olabel'), primitive.copy().arcsort('ilabel'))
    else:
        f = fstlib.compose(primitive.copy().arcsort('olabel'), words.arcsort('ilabel'))
    return fstlib.determinize(f.project(side).rmepsilon()).minimize()


def _total(f):
    """Shortest distance from the start state to any final state (inf if none)."""
    if f.start() < 0:
        return math.inf
    d = fstlib.shortestdistance(f, reverse=True)
    return float(d[f.start()]) if f.start() < len(d) else math.inf


def _words(f, sym):
    """Every word of an acyclic acceptor, with its path cost."""
    out = {}
    stack = [(f.start(), '', 0.0)]
    while stack:
        state, word, cost = stack.pop()
        final = float(f.final(state))
        if math.isfinite(final):
            out[word] = min(out.get(word, math.inf), cost + final)
        for arc in f.arcs(state):
            stack.append((arc.nextstate, word + sym.find(arc.ilabel), cost + float(arc.weight)))
    return out
