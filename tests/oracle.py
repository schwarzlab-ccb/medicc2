import math
import numpy as np
import medicc

def get_parts():
    sym = medicc.create_symbol_table(8, 'X')
    p = medicc.create_copynumber_fst(symbol_table=sym, sep='X', enable_wgd=True,
                                  max_num_wgds=3, wgd_x2=True, output_all=True)
    p['W1'] = medicc.factory.create_1step_WGD_fst(sym, 'X', wgd_cost=1, minimize=True, wgd_x2=True)
    # Keep L_LOH and G separate. This independently checks the detmin NOWGD path.
    return p, sym

class Primitive:
    def __init__(self, wrapper):
        f = wrapper.fst
        self.start = f.start()
        self.final = {s: float(f.final(s)) for s in f.states() if math.isfinite(float(f.final(s)))}
        self.transitions = {}
        for s in f.states():
            for a in f.arcs(s):
                assert a.ilabel and a.olabel
                self.transitions.setdefault((s, a.ilabel, a.olabel), []).append((a.nextstate, float(a.weight)))

    def score(self, a, b):
        states = {self.start: 0.0}
        for x, y in zip(a, b):
            nxt = {}
            for s, v in states.items():
                for t, w in self.transitions.get((s, x, y), ()):
                    nv = v + w
                    if nv < nxt.get(t, math.inf): nxt[t] = nv
            states = nxt
        return min((v + self.final.get(s, math.inf) for s, v in states.items()), default=math.inf)

def minplus(a, b):
    assert a.shape[1] == b.shape[0]
    out = np.full((a.shape[0], b.shape[1]), np.inf)
    for k in range(a.shape[1]):
        np.minimum(out, a[:, k:k+1] + b[k:k+1, :], out=out)
    return out

def relation(primitive, strings, sym):
    labels = [[sym.find(c) for c in s] for s in strings]
    return np.array([[primitive.score(a, b) for b in labels] for a in labels])

def optimum(ss,T,a,b):
    ix={str(s):i for i,s in enumerate(ss)}
    av=np.full(len(ss),np.inf);bv=av.copy()
    for s,w in a.items():av[ix[s]]=w
    for s,w in b.items():bv[ix[s]]=w
    costs=np.min(T+av[None,:],axis=1)+np.min(T+bv[None,:],axis=1)
    low=np.min(costs)
    return {} if not np.isfinite(low) else {str(s):int(w) for s,w in zip(ss,costs) if w==low}
