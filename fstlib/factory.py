'''
Created on 21 Jan 2010

@author: rfs
'''


import numpy as np
import collections
from itertools import combinations
import fstlib

""" Creates special transducers for special occasions ;) """
    
EDNAFULL = {('a','a'):5,
            ('t','t'):5,
            ('g','g'):5,
            ('c','c'):5,
            ('a','t'):-4,
            ('a','g'):-4,
            ('a','c'):-4,
            ('c','a'):-4,
            ('c','g'):-4,
            ('c','t'):-4,
            ('g','a'):-4,
            ('g','c'):-4,
            ('g','t'):-4,
            ('t','a'):-4,
            ('t','g'):-4,
            ('t','c'):-4}

EDITDNA = {('a','a'):0,
            ('t','t'):0,
            ('g','g'):0,
            ('c','c'):0,
            ('a','t'):1,
            ('a','g'):1,
            ('a','c'):1,
            ('c','a'):1,
            ('c','g'):1,
            ('c','t'):1,
            ('g','a'):1,
            ('g','c'):1,
            ('g','t'):1,
            ('t','a'):1,
            ('t','g'):1,
            ('t','c'):1}

#DNAEXP = {('A', 'A'): 2.8512358587881401, ('A', 'T'): 3.0396830090087774, ('T', 'T'): 1.8974875643790847, ('T', 'A'): 3.0396830090087774, ('C', 'A'): 3.7190999967568206, ('C', 'T'): 3.2662887680514583, ('G', 'G'): 2.8028092648826655, ('T', 'C'): 3.2662887680514583, ('G', 'A'): 3.6155593178159804, ('G', 'T'): 2.8398505365630147, ('A', 'G'): 3.6155593178159804, ('C', 'G'): 3.3569853295222702, ('C', 'C'): 2.4181654845880867, ('T', 'G'): 2.8398505365630147, ('G', 'C'): 3.3569853295222702, ('A', 'C'): 3.7190999967568206}
DNAEXP = {('A', 'A'): 2.3797654605525818, ('A', 'T'): 3.1532202470020141, ('T', 'T'): 2.3500179409534674, ('T', 'A'): 3.1532202470020141, ('C', 'A'): 3.1664024667929485, ('C', 'T'): 3.1572330934783825, ('G', 'G'): 2.392298397615908, ('T', 'C'): 3.1572330934783825, ('G', 'A'): 3.1828086852453756, ('G', 'T'): 3.1610199086227939, ('A', 'G'): 3.1828086852453756, ('C', 'G'): 3.1714878096208929, ('C', 'C'): 2.4101012294965733, ('T', 'G'): 3.1610199086227939, ('G', 'C'): 3.1714878096208929, ('A', 'C'): 3.1664024667929485}

def _adjust_scoring_matrix(scoring_matrix, invert = True):
    """the BioPython scoring matrices are not mirrored, i.e. if A->T is defined
    T->A isn't. This function fixes that. It additionally also multiplies all scores
    by -1 if 'invert = True' and turns all into lower case letters."""
    #unique_keys = set(map(lambda x:tuple(sorted(x)),scoring_matrix.keys()))
    new_scoring_matrix={}
        
    for key in scoring_matrix.keys():
        new_key = tuple([k.lower() for k in key])
        new_scoring_matrix[new_key] = scoring_matrix[key] * (-1 if invert else 1)
        new_scoring_matrix[new_key[::-1]] = new_scoring_matrix[new_key]
    
    return new_scoring_matrix
    
def create_three_state_global_alignment_fst(scoring_matrix, gap_open, gap_extend, invert_matrix = False, gap_symbol='-'):
    """creates a transducer for global alignment with
    specified gap open, extend and scoring costs. This FST 
    has 3 states and distinguishes between inserts and deletes """

    scoring_matrix = _adjust_scoring_matrix(scoring_matrix, invert_matrix)

    alphabet = set([i[0] for i in scoring_matrix.keys()])
    symbol_table = fstlib.SymbolTable()
    symbol_table.add_symbol(gap_symbol)
    for s in alphabet:
        symbol_table.add_symbol(s)
        
    myfst = fstlib.Fst(arc_type='standard')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)
    myfst.add_states(3)
    myfst.set_start(0)
    myfst.set_final(0,0)
    myfst.set_final(1,0)
    myfst.set_final(2,0)

    # starting state is match state
    # create match transitions to self
    myfst.add_arcs(0, [(s, t, scoring_matrix[(s,t)], 0) for s in alphabet for t in alphabet])

    # create gap open transitions to gap state [delete]
    myfst.add_arcs(0, [(s, gap_symbol, gap_open, 1) for s in alphabet])        
        
    # create gap open transitions to gap state [insert]
    myfst.add_arcs(0, [(gap_symbol, s, gap_open, 2) for s in alphabet])
        
    # create self transitions to gap states (gap extend) [delete]
    myfst.add_arcs(1, [(s, gap_symbol, gap_extend, 1) for s in alphabet])        
        
    # create self transitions to gap states (gap extend) [insert]
    myfst.add_arcs(2, [(gap_symbol, s, gap_extend, 2) for s in alphabet])

    # create match transitions back to match state
    myfst.add_arcs(2, [(s,t, scoring_matrix[(s,t)], 0) for s in alphabet for t in alphabet])
    myfst.add_arcs(1, [(s,t, scoring_matrix[(s,t)], 0) for s in alphabet for t in alphabet])
        
    # create delete to delete transitions
    myfst.add_arcs(1, [(gap_symbol, s, gap_open, 2) for s in alphabet])
    myfst.add_arcs(2, [(s, gap_symbol, gap_open, 1) for s in alphabet])
        
    return myfst, symbol_table

def create_kgram_fst(symbol_table, k, gap_cost=None):
    st = dict(symbol_table)
    myfst = fstlib.Fst(arc_type='log')
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)
    myfst.add_states(k+1)
    myfst.set_start(0)
    myfst.set_final(k,0)
    gap = st[0]
    
    # state 0
    myfst.add_arcs(0, [(s, gap, 0, 0) for key,s in symbol_table if key!=0])
        
    for i in range(k):
        myfst.add_arcs(i, [(s,s, 0, i+1) for key,s in symbol_table if key!=0])
        if gap_cost is not None and i>0: ## gappy k-gram kernel
            myfst.add_arcs(i, [(s, gap, -np.log(gap_cost), i) for key,s in symbol_table if key!=0])

    # final state
    myfst.add_arcs(k, [(s, gap, 0, k) for key,s in symbol_table if key!=0])

    return myfst

def from_string(seq, final_weight=None, arc_weight=None, arc_type = 'standard', isymbols=None, osymbols=None):
    myfst = fstlib.Fst(arc_type=arc_type)
    startid = myfst.add_state()
    myfst.set_start(startid)
    myfst.set_input_symbols(isymbols)
    myfst.set_output_symbols(osymbols)
    for i,s in enumerate(seq):
        myfst.add_state()
        if arc_weight is None:
            arc_weight = fstlib.Weight.one(myfst.weight_type())
        else:
            if not isinstance(arc_weight, fstlib.Weight):
                arc_weight = fstlib.Weight(myfst.weight_type(), arc_weight)
        myfst.add_arc(i, (s, s, arc_weight, i+1))
    if final_weight is None:
        final_weight = fstlib.Weight.one(myfst.weight_type())
    else:
        if isinstance(final_weight, fstlib.Weight):
            pass
        else:
            final_weight = fstlib.Weight(myfst.weight_type(), final_weight)
    myfst.set_final(i+1, final_weight)
    return myfst

def from_array(seq, final_weight=None, arc_weight=None, arc_type = 'standard', symbols=None):
    myfst = fstlib.Fst(arc_type=arc_type)
    startid = myfst.add_state()
    myfst.set_start(startid)
    myfst.set_input_symbols(symbols)
    myfst.set_output_symbols(symbols)
    #for i,s in enumerate(seq):
    last_id = startid
    for i in range(seq.shape[1]): ## for each column
        current_id = myfst.add_state()        
        for j in range(seq.shape[0]): ## for each row
            entry = seq[j, i]
            if seq.ndim == 3:
                s, t = entry
            elif seq.ndim == 2:
                s = t = entry
            else:
                raise FSTFactoryError("Wrong number of dimensions in input array")
            if arc_weight is None:
                arc_weight = fstlib.Weight.one(myfst.weight_type())
            else:
                if not isinstance(arc_weight, fstlib.Weight):
                    arc_weight = fstlib.Weight(myfst.weight_type(), arc_weight)
            myfst.add_arc(last_id, (s, t, arc_weight, current_id))
        last_id = current_id

    if final_weight is None:
        final_weight = fstlib.Weight.one(myfst.weight_type())
    else:
        if isinstance(final_weight, fstlib.Weight):
            pass
        else:
            final_weight = fstlib.Weight(myfst.weight_type(), final_weight)
    myfst.set_final(last_id, final_weight)
    return myfst

def from_count_matrix(count_matrix, symbol_table, keys=None, arc_type = 'log', normalize=True):
    """ Creates a simple linear prob. FSA from a matrix of counts.
    The count matrix has a number of rows equal to the number of symbols and a number of
    columns equal to the length of the sequence. The keys parameter maps row indices in the count matrix
    to symbol table keys, i.e. keys[i] defines the symbol table key of the i-th row in the count matrix.
    If normalize is given, all column counts are normalized to sum to 1 (defines probability distribution). """
    st = dict(symbol_table)
    if keys is None:
        keys = list(st.keys())
    if count_matrix.shape[0] != len(keys):
        raise FSTFactoryError("Number of rows of count_matrix does not match symbol table or keys.")
    
    seq_length = count_matrix.shape[1]

    myfst = fstlib.Fst(arc_type)
    myfst.set_input_symbols(symbol_table)
    myfst.set_output_symbols(symbol_table)
    myfst.add_states(seq_length + 1)	
    myfst.set_start(0)
    myfst.set_final(seq_length, 0)

    count_sum = count_matrix.sum(0)
    if np.any(count_sum==0):
        raise FSTFactoryError("Column sumns must not be zero.")
                
    for i in range(0, seq_length):
        which = np.nonzero(count_matrix[:,i]>0)[0]
        for j in which:
            key = keys[j]
            prob = count_matrix[j,i]
            if normalize:
                prob = float(prob) / count_sum[i]
            if prob > 0:
                prob = - np.log(prob)
                myfst.add_arc(i, (key, key, prob, i+1))
        
    return myfst

def create_fixed_length_fsa_from_symbol_table(symbol_table, length, epsilons=False, arc_type='standard'):
    n = length+1
    fst = fstlib.Fst(arc_type=arc_type)
    fst.add_states(n)
    fst.set_start(0)
    fst.set_final(n-1, 0)
    fst.set_input_symbols(symbol_table)
    fst.set_output_symbols(symbol_table)
    for i in range(n-1):
        if epsilons:
            fst.add_arcs(i, [(s,s,0,i+1) for _,s in symbol_table])
        else:
            fst.add_arcs(i, [(s,s,0,i+1) for k,s in symbol_table if k != 0])

    return fst

class FSTFactoryError(Exception):
    def __init__(self, msg):
        self.message = msg

    def __str__(self):
        return self.message

