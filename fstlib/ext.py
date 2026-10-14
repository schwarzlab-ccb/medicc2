import logging

import pandas as pd
import numpy as np

import fstlib
import fstlib.cext.ops
from fstlib.cext import pywrapfst

logger = logging.getLogger(__name__)

try:
    from fstlib.cext.ops import *
    logger.info('Using accelerated C operations.')
except ImportError:
    logger.warn('Accelerated C operations not available.')

def align(model, ifst1, ifst2):
    ofst = fstlib.cext.ops.align(model.fst, ifst1.fst, ifst2.fst)
    return fstlib.Fst(ofst)

def decompose(ifsa, ofsa, event_fst, n_events, return_all=False):
    C = [None, ] * n_events ## cliques

    ## initial beliefs
    C[0] = ifsa * event_fst
    for i in range(1, n_events-1):
        C[i] = event_fst
    C[n_events-1] = event_fst * ofsa

    ## messages
    fwd = [None, ] * (n_events-1)
    bwd = [None, ] * (n_events-1)

    fwd[0] = fstlib.determinize(fstlib.project(C[0], 'output')).minimize()
    for i in range(1, n_events-1):
        fwd[i] = fstlib.determinize(fstlib.project(fwd[i-1] * C[i], 'output')).minimize()
    
    bwd[n_events-2] = fstlib.determinize(fstlib.project(C[n_events-1], 'input')).minimize()
    for i in range(n_events-3, -1, -1):
        bwd[i] = fstlib.determinize((C[i+1] * bwd[i+1]).project('input')).minimize()

    ## final beliefs
    B = [None, ] * n_events
    B[0] = C[0] * bwd[0]
    for i in range(1, n_events-1):
        B[i] = fwd[i-1] * C[i] * bwd[i]
    B[n_events-1] = fwd[n_events-2] * C[n_events-1]

    if return_all:
        return C, fwd, bwd, B
    else: ## only return final beliefs
	    return B

def encode_determinize_minimize(ifst, delta=1e-6, minimize=True):
    em = fstlib.EncodeMapper(arc_type=ifst.arc_type(), encode_labels=True, encode_weights=False)
    ofst = ifst.copy().encode(em)
    ofst = fstlib.determinize(ofst, delta=delta)
    if minimize:
        ofst.minimize()
    ofst.decode(em)
    return ofst

def info(ifst, name=None):
    """ Outputs some info stats """
    info = {}
    info['fst_type'] = ifst.fst_type()
    info['arc_type'] = ifst.arc_type()
    info['nstates'] = ifst.num_states()
    info['nfinalstates'] = np.sum([ifst.is_final(s) for s in ifst.states()])
    info['narcs'] = 0
    info['ninputeps'] = 0
    info['noutputeps'] = 0
    for state in ifst.states():
        info['narcs'] += ifst.num_arcs(state)
        info['ninputeps'] += ifst.num_input_epsilons(state)
        info['noutputeps'] += ifst.num_output_epsilons(state)
    info['acceptor'] = ifst.properties(fstlib.FstProperties.ACCEPTOR, True) == fstlib.FstProperties.ACCEPTOR
    info['input_deterministic'] = ifst.properties(fstlib.FstProperties.I_DETERMINISTIC, True) == fstlib.FstProperties.I_DETERMINISTIC
    info['output_deterministic'] = ifst.properties(fstlib.FstProperties.O_DETERMINISTIC, True) == fstlib.FstProperties.O_DETERMINISTIC
    info['input_label_sorted'] = ifst.properties(fstlib.FstProperties.I_LABEL_SORTED, True) == fstlib.FstProperties.I_LABEL_SORTED
    info['output_label_sorted'] = ifst.properties(fstlib.FstProperties.O_LABEL_SORTED, True) == fstlib.FstProperties.O_LABEL_SORTED
    info['cyclic'] = ifst.properties(fstlib.FstProperties.CYCLIC, True) == fstlib.FstProperties.CYCLIC
    info['topsorted'] = ifst.properties(fstlib.FstProperties.TOP_SORTED, True) == fstlib.FstProperties.TOP_SORTED
    info['accessible'] = ifst.properties(fstlib.FstProperties.ACCESSIBLE, True) == fstlib.FstProperties.ACCESSIBLE
    info['coaccessible'] = ifst.properties(fstlib.FstProperties.COACCESSIBLE, True) == fstlib.FstProperties.COACCESSIBLE
    info['weighted'] = ifst.properties(fstlib.FstProperties.WEIGHTED, True) == fstlib.FstProperties.WEIGHTED

    df = pd.DataFrame.from_dict(info, orient='index', columns=[name if name is not None else id(ifst)])
    return df

def kernel_score(model, ifst1, ifst2):
    if model.arc_type()=='standard':
        distance = fstlib.cext.ops.kernel_score_std(model.fst, ifst1.fst, ifst2.fst)
    elif model.arc_type()=='log':
        distance = fstlib.cext.ops.kernel_score_log(model.fst, ifst1.fst, ifst2.fst)
    else:
        raise FSTlibExtError('Kernel score not implemented for %s semiring' % model.arc_type())
    return distance

def multi_score(fst1, fst2, fst3, ifst1, ifst2):
    if fst1.arc_type()=='standard':
        distance = fstlib.cext.ops.multi_score_std(fst1.fst, fst2.fst, fst3.fst, ifst1.fst, ifst2.fst)
    else:
        raise FSTlibExtError('Multi score not implemented for %s semiring' % fst1.arc_type())
    return distance

def multi_kernel_score(fst1, fst2, fst3, fst4, ifst1, ifst2):
    if fst1.arc_type()=='standard':
        distance = fstlib.cext.ops.multi_kernel_score_std(fst1.fst, fst2.fst, fst3.fst, fst4.fst, ifst1.fst, ifst2.fst)
    else:
        raise FSTlibExtError('Multi kernel score not implemented for %s semiring' % fst1.arc_type())
    return distance

def normalize(ifst, inplace=False):
    """Normalizes fst so that all outgoing transitions sum to 1."""
    
    ## convert to real if necessary
    if inplace:
        ofst = ifst
    else:
        ofst = ifst.copy()
        
    ofst.weight_map(fstlib.algos.map_log_to_real)

    for state in ofst.states():
        total_arcweights = np.sum([float(arc.weight) for arc in ofst.arcs(state)])
        mai = ofst.mutable_arcs(state)
        for arc in mai:
            arc.weight = fstlib.Weight(ofst.weight_type(), float(arc.weight) / total_arcweights)
            mai.set_value(arc)

    ## convert back
    ofst.weight_map(fstlib.algos.map_real_to_log)
    return ofst

def normalize_alphabet(ifst, inplace=False):
    """Normalizes fst so that outgoing transition weights of the same input symbol sum to 1"""
    
    ## convert to real if necessary
    if inplace:
        ofst = ifst
    else:
        ofst = ifst.copy()
        
    ofst.weight_map(fstlib.tools.neglog_to_real)

    for state in ofst.states():
        arcweights = [(arc.ilabel, float(arc.weight)) for arc in ofst.arcs(state)]
        labelweights = pd.DataFrame.from_records(arcweights, columns=['ilabel','weight']).groupby('ilabel').sum()
        mai = ofst.mutable_arcs(state)
        for arc in mai:
            labelsum = labelweights.weight.loc[arc.ilabel]
            arc.weight = fstlib.Weight(ofst.weight_type(), float(arc.weight) / labelsum)
            mai.set_value(arc)

    ## convert back
    ofst.weight_map(fstlib.tools.real_to_neglog)
    return ofst

def project(ifst, project_type='input'):
    """ Provides a non-destructive project operation, which pywrapfst doesn't have. """
    ofst = ifst.copy()
    ofst.project(project_type)
    return ofst

def read(source):
    newfst = pywrapfst.Fst.read(source)
    return fstlib.Fst(newfst)

def score(model, ifst1, ifst2):
    if model.arc_type()=='standard':
        distance = fstlib.cext.ops.score_std(model.fst, ifst1.fst, ifst2.fst)
    elif model.arc_type()=='log':
        distance = fstlib.cext.ops.score_log(model.fst, ifst1.fst, ifst2.fst)
    else:
        pass
    return distance

def weight_map(ifst, func):
    newfst = ifst.copy()
    newfst.weight_map(func)
    return newfst


class FSTlibExtError(Exception):
    pass