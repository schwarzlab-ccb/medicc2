import logging

import fstlib
import numpy as np

from medicc import tools

SYMBOL_GAP = 'gap'

# prepare logger
logger = logging.getLogger(__name__)


def create_symbol_table(max_cn=8, separator='X'):
    alphabet = [tools.int2hex(x) for x in range(max_cn+1)]
    symbol_table = fstlib.SymbolTable()
    symbol_table.add_symbol(SYMBOL_GAP)
    for s in alphabet:
        symbol_table.add_symbol(s)

    if separator is not None and separator != '':
        symbol_table.add_symbol(separator)
    return symbol_table


def _get_int_cns_from_symbol_table(symbol_table, separator='X'):
    cns = {x[1]: tools.hex2int(x[1]) for x in symbol_table if x[0] != 0 and x[1] != separator}
    return cns


def create_1step_del_fst(symbol_table, separator='X', exclude_zero=False, w_stay=0, w_open=1, w_extend=0):

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)

    myfst.add_states(2)
    myfst.set_start(0)
    myfst.set_final(0, 0)
    myfst.set_final(1, 0)

    myfst.add_arcs(0, [(s, s, w_stay, 0) for s in cns.keys()])
    if separator is not None and separator != '':
        myfst.add_arc(0, (separator, separator, 0, 0))  # add separator
    if exclude_zero:
        myfst.add_arcs(0, [(s, t, w_open, 1) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == 1 and t != '0'])  # 0->1
        myfst.add_arcs(1, [(s, t, w_extend, 1) for s in cns.keys() for t in cns.keys()
                           if (cns[s]-cns[t]) == 1 and t != '0'])  # extend an open window
    else:
        myfst.add_arcs(0, [(s, t, w_open, 1) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == 1])  # 0->1
        myfst.add_arcs(1, [(s, t, w_extend, 1) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == 1])  # extend an open window
    if '0' in cns.keys():
        myfst.add_arc(1, ('0', '0', 0, 1))
    myfst.add_arcs(1, [(s, s, w_stay, 0) for s in cns.keys() if s != '0'])  # return
    if separator is not None and separator != '':
        myfst.add_arc(1, (separator, separator, 0, 0))

    return myfst


def create_loh_fst(symbol_table, separator='X'):

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)

    myfst.add_states(9)
    myfst.set_start(0)
    for i in range(9):
        myfst.set_final(i, 0)

    ## match 0 -> 0
    myfst.add_arcs(0, [(s, s, 0, 0) for s in cns.keys()])
    if separator is not None and separator != '':
        myfst.add_arc(0, (separator, separator, 0, 0))  # add separator

    ## others
    for state in range(1, 9):
        cost = state
        myfst.add_arcs(0, [(s, t, cost, state) for s in cns.keys()
                           for t in cns.keys() if (cns[s]-cns[t]) == cost and t == '0'])
        myfst.add_arcs(state, [(s, t, 0, state) for s in cns.keys()
                               for t in cns.keys() if (cns[s]-cns[t]) <= cost and t == '0'])
        myfst.add_arcs(state, [(s, t, cns[s]-cns[t]-cost, cns[s]-cns[t]) for s in cns.keys()
                               for t in cns.keys() if (cns[s]-cns[t]) > cost and t == '0'])
        myfst.add_arcs(state, [(s, s, 0, 0) for s in cns.keys()])
        if separator is not None and separator != '':
            myfst.add_arc(state, (separator, separator, 0, 0))

    myfst = fstlib.encode_determinize_minimize(myfst)
    
    return myfst


def create_1step_WGD_fst(symbol_table, separator='X', wgd_cost=1, minimize=True, wgd_x2=False,
                         total_cn=False):

    if total_cn:
        wgd_distance = 2.
    else:
        wgd_distance = 1.

    cns = _get_int_cns_from_symbol_table(symbol_table, separator)

    W = fstlib.Fst()
    W.set_input_symbols(symbol_table)
    W.set_output_symbols(symbol_table)
    W.add_states(3)
    W.set_start(0)
    W.set_final(0, 0)
    W.set_final(1, 0)
    W.set_final(2, 0)
    if wgd_x2:
        W.add_arcs(0, [(s, t, wgd_cost, 1) for s in cns.keys()
                    for t in cns.keys() if (s != '0') and ((cns[t]/cns[s]) == 2.)])
    else:
        W.add_arcs(0, [(s, t, wgd_cost, 1) for s in cns.keys()
                    for t in cns.keys() if (s != '0') and ((cns[t]-cns[s]) == wgd_distance)])

    W.add_arc(1, ('0', '0', 0, 1))
    if wgd_x2:
        W.add_arcs(1, [(s, t, 0, 1) for s in cns.keys()
                    for t in cns.keys() if (s != '0') and ((cns[t]/cns[s]) == 2.)])
    else:
        W.add_arcs(1, [(s, t, 0, 1) for s in cns.keys()
                    for t in cns.keys() if (s != '0') and ((cns[t]-cns[s]) == wgd_distance)])
    W.add_arc(0, ('0', '0', 0, 0))
    if separator is not None and separator != '':
        W.add_arc(0, (separator, separator, 0, 0))
        W.add_arc(1, (separator, separator, 0, 1))
    W.add_arc(1, ('0', '0', 0, 1))

    W.add_arcs(0, [(s, s, 0, 2) for s in cns.keys() if s != '0'])
    W.add_arcs(2, [(s, s, 0, 2) for s in cns.keys()])
    W.add_arc(2, ('0', '0', 0, 2))
    if separator is not None and separator != '':
        W.add_arc(2, (separator, separator, 0, 2))
    W.arcsort('olabel')
    if minimize:
        W = fstlib.encode_determinize_minimize(W)

    return W


def create_nstep_fst(n, one_step_fst, minimize=True):
    # Extend 1step FST
    nstep_fst = one_step_fst
    nstep_fst.arcsort(sort_type='olabel')

    for _ in range(n):
        nstep_fst = fstlib.compose(nstep_fst, one_step_fst)
        if minimize:
            nstep_fst = fstlib.encode_determinize_minimize(nstep_fst)
        nstep_fst.arcsort(sort_type='olabel')

    return nstep_fst


def create_copynumber_fst(symbol_table, sep='X', enable_wgd=False, wgd_cost=1, 
                          max_num_wgds=3, wgd_x2=False, output_all=False, total_cn=False,
                          exact_nowgd=True, max_pre_wgd_losses=8, exact_wgd=False):
    """ Creates the tree FST T which computes the asymmetric MED.
    The current creation is based on a trade-off. In the absence of WGDs (and the default "exact_nowgd=True"),
    the FST is exact wr.t. combined LOH-losses (i.e 21 -> 10 is counted as one event), however,
    in the presence of WGDs, for performance reasons losses and LOHs are counted separately.

    For `exact_wgd=True`, the FST is exact w.r.t. combined LOH-losses even in the presence of WGDs 
    but due to the exponential growth of the FST, it is not feasible in most cases. (Example where 
    this is required: 22X1X1X1 -> 10X2X2X2, 4 events with exact_wgd and 5 otherwise.)

    For `exact_nowgd=False` the legacy version is created which never combines losses and LOHs.
    """
    n = len(_get_int_cns_from_symbol_table(symbol_table, sep))

    L_1step = create_1step_del_fst(symbol_table, sep, exclude_zero=True)
    L = create_nstep_fst(n-1, L_1step)
    LG = fstlib.encode_determinize_minimize(L*~L)
    G = ~L
    L_LOH_1step = create_1step_del_fst(symbol_table, sep, exclude_zero=False)
    L_LOH = create_nstep_fst(max_pre_wgd_losses-1, L_LOH_1step)
    LOH = create_loh_fst(symbol_table, sep)
    
    if enable_wgd:
        W1step = create_1step_WGD_fst(symbol_table, sep, wgd_cost=wgd_cost,
                                      minimize=False, wgd_x2=wgd_x2, total_cn=total_cn)
        n = int(min(max_num_wgds, np.floor(np.log2(n))))
        if n > 1:
            W = create_nstep_fst(n-1, W1step)
        else:
            W = W1step
        if exact_wgd:
            T = L_LOH * W * LG
        elif exact_nowgd:
            T = ((LOH * W * LG) + fstlib.encode_determinize_minimize(L_LOH * G)).rmepsilon()
        else:
            # legacy version
            T = LOH * W * LG
    else:
        if exact_nowgd:
            T = fstlib.encode_determinize_minimize(L_LOH*G)
        else:
            # legacy version
            T = LOH * LG
        W = None

    if output_all:
        return {'T': T, 'LOH': LOH, 'W': W, 'L': L, 'L_LOH': L_LOH, 'G': G, 'LG': LG}
    else:
        return T
