"""Literal-FST test oracle, built with fstlib (OpenFST) from MEDICC's factory primitives.

It builds the deployed edge model literally:

    P = A + ADA + ADADA,   A = L_LOH o G,   D: c -> min(2c, 8), X -> X (mandatory, cost 1)

i.e. one chain A (D A)^j per WGD count j = 0, 1, 2, and T(x, y) is the cheapest path
through any of the chains. Each WGD costs 1 once per genome (the native engine adds k + j
for a branch with k and j WGDs). OpenFST is the reference; nothing here shares code with
the native engines.

The chains are never composed into one transducer (ADADA is far too large).
Instead a weighted word language is pushed through one primitive at a time, and
determinized after each step, so only the words reachable from the input remain.
The results of the chains are then united; determinizing in the tropical semiring keeps
the minimum cost of each word, which is the "+" of P.

fstlib operators: `a * b` composes a with b (it arc-sorts a by output label first); on two
acceptors that is their intersection. `a + b` is the union (of a copy of a).
"""
from collections.abc import Mapping, Sequence
from functools import lru_cache
import math

import fstlib
import medicc

MAX_CN = 8
SEPARATOR = 'X'
MAX_WGDS = 2

Chain = Sequence[fstlib.Fst]


def get_parts():
    """Return (primitives, symbol table): MEDICC's L_LOH and G, plus the doubling D."""
    sym = medicc.create_symbol_table(MAX_CN, SEPARATOR)
    p = medicc.create_copynumber_fst(symbol_table=sym, sep=SEPARATOR, enable_wgd=False,
                                     output_all=True)
    p = {'L_LOH': p['L_LOH'], 'G': p['G'], 'D': doubling(sym, cost=1)}
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


def chain(parts, wgds):
    """A (D A)^wgds as a list of primitives: L_LOH, G, D, L_LOH, G, ..."""
    a = [parts['L_LOH'], parts['G']]
    return a + [parts['D'], *a] * wgds


class Oracle:
    """Exact scores and ancestor languages for short words, computed with OpenFST."""

    def __init__(self, wgds=MAX_WGDS):
        self.parts, self.sym = get_parts()
        self.model = [chain(self.parts, j) for j in range(wgds + 1)]  # A, ADA, ADADA

    def language(self, words: Mapping[str, float]) -> fstlib.Fst:
        """Weighted acceptor: each word with its cost."""
        return _union([fstlib.factory.from_string(word, final_weight=cost, isymbols=self.sym,
                                                  osymbols=self.sym)
                       for word, cost in words.items()])

    def children(self, parents: fstlib.Fst, model: Sequence[Chain] = None) -> fstlib.Fst:
        """Every child word y with cost min_x parents[x] + T(x, y)."""
        out = []
        for c in model or self.model:
            words = parents
            for p in c:
                words = _step(words, p, 'output')
            out.append(words)
        return _union(out)

    def ancestors(self, children: fstlib.Fst, model: Sequence[Chain] = None) -> fstlib.Fst:
        """Every ancestor word z with cost min_y T(z, y) + children[y]."""
        out = []
        for c in model or self.model:
            words = children
            for p in reversed(c):
                words = _step(words, p, 'input')
            out.append(words)
        return _union(out)

    def score(self, x: str, y: str, model: Sequence[Chain] = None) -> float:
        """Cheapest path from x to y through `model` (default: the edge model P)."""
        return _total(self.children(self.language({x: 0}), model) * self.language({y: 0}))

    @lru_cache(maxsize=None)
    def directed(self, x: str, y: str) -> float:
        """T(x, y): cheapest edge history from parent x to child y."""
        return self.score(x, y)

    @lru_cache(maxsize=None)
    def kernel(self, x: str, y: str) -> float:
        """K(x, y) = min_z T(z, x) + T(z, y): the pair distance."""
        return self.pair(x, y)

    def pair(self, x: str, y: str, left: Sequence[Chain] = None, right: Sequence[Chain] = None) -> float:
        """min_z left(z, x) + right(z, y). Both default to the edge model P."""
        zx = self.ancestors(self.language({x: 0}), left)
        zy = self.ancestors(self.language({y: 0}), right)
        return _total(zx * zy)

    def optimum(self, a: Mapping[str, float], b: Mapping[str, float]) -> dict[str, int]:
        """All optimal ancestors z of two weighted child languages, with their costs.

        cost(z) = min_w (T(z, w) + a[w]) + min_w (T(z, w) + b[w]); only the minimum-cost z are kept.
        """
        joint = self.ancestors(self.language(a)) * self.ancestors(self.language(b))
        if not math.isfinite(_total(joint)):
            return {}
        best = fstlib.determinize(joint)
        best.prune(weight=0)  # keep every word that ties the optimum
        return {word: int(cost) for word, cost in _words(best, self.sym).items()}


def _union(fsts):
    """One acceptor of all words of `fsts`; a word in several keeps its minimum cost."""
    out = None
    for f in fsts:
        out = f if out is None else out + f
    return fstlib.determinize(out.rmepsilon()).minimize()


def _step(words, primitive, side):
    """Push a weighted word language through one primitive.

    side='output': words are inputs of `primitive`; keep its outputs (parent -> child).
    side='input':  words are outputs of `primitive`; keep its inputs (child -> parent).
    """
    f = words * primitive if side == 'output' else primitive * words
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
