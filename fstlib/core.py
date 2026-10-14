""" Python class for finite-state transducers over the real or log semiring"""

import logging
import math
import os
import sys
import tempfile
from io import StringIO
import collections

import numpy as np
import pandas as pd

import fstlib
import fstlib.algos
from fstlib.cext import pywrapfst

logger = logging.getLogger(__name__)

class Fst:
    """ Class that describes a finite-state transducer.
    """
    
    def __init__(self, arc_type = fstlib.Semiring.TROPICAL):
        if isinstance(arc_type, pywrapfst.MutableFst):
            self.fst = arc_type
        elif isinstance(arc_type, Fst):
            self.fst = arc_type.fst
        else:
            self.fst = pywrapfst.VectorFst(arc_type=arc_type)

    def __repr__(self):
        return self.fst.__repr__()

    def _repr_png_(self):
        if self.num_states() > 50: ## don't draw, too large, too slow
            logger.info('Too many states to draw')
            return
        g = self.to_graphviz(width=12, height=14)
        g.format='png'
        return g.pipe()

    def _repr_svg_(self):
        if self.num_states() > 50: ## don't draw, too large, too slow
            logger.info('Too many states to draw')
            return
        g = self.to_graphviz(width=12, height=14)
        g.format='svg'
        return str(g.pipe(), 'utf-8')

    def __str__(self):
        return self.fst.__str__()

    def __mul__(self, other):
        return fstlib.compose(self.arcsort('olabel'), other)

    def __add__(self, other):
        return self.copy().union(other)

    def __invert__(self):
        return self.copy().invert()
    
    def add_arc(self, state, arc):
        """ Arc can be pywrapfst.Arc or 4-tuple (ilabel, olabel, weight, nextstate).
         An extended version of pywrapfst.MutableFst.add_arc() that can handle tuples as input."""
        if isinstance(arc, pywrapfst.Arc):
            pass
        else: ## assume tuple / unpackable
            ### adds an arc object expecting the right weight type and trying
            ### to resolve symbolic (non-integer) labels via the symbol table.
            ilabel, olabel, weight, nextstate = arc
            isyms = self.input_symbols()
            osyms = self.output_symbols()
            if isinstance(ilabel, str):
                if isyms is not None:
                    isym = isyms.find(ilabel)
                else:
                    raise FSTError('ilabel is String but no input symbol table available.')
            else:
                isym = ilabel

            if isinstance(olabel, str):
                if osyms is not None:
                    osym = osyms.find(olabel)
                else:
                    raise FSTError('olabel is String but no output symbol table available.')
            else:
                osym = olabel

            if isinstance(weight, pywrapfst.Weight):
                pass
            else:
                weight = pywrapfst.Weight(self.weight_type(), weight)
            arc = pywrapfst.Arc(isym, osym, weight, nextstate)
        
        self.fst.add_arc(state, arc)
        return self

    def add_arcs(self, state, iterable_of_arcs):
        """ Add an iterable of arcs or 4-tuples to a state.
        An extension to pywrapfst.MutableFst."""
        for arc in iterable_of_arcs:
            self.add_arc(state, arc)
        return self

    def add_state(self):
        return self.fst.add_state()

    def add_states(self, n):
        self.reserve_states(n)
        ids = np.array([self.add_state() for i in range(n)], dtype=np.int32)
        return ids
        
    def arc_type(self):
        return self.fst.arc_type()

    def arcs(self, state):
        return self.fst.arcs(state)
    
    def arcsort(self, sort_type='ilabel'):
        self.fst.arcsort(sort_type)
        return self

    def closure(self, closure_type = "star"):
        self.fst.closure(closure_type)
        return self

    def concat(self, ifst):
        self.fst.concat(ifst.fst)
        return self

    def connect(self):
        self.fst.connect()
        return self

    def copy(self):
        newfst = self.fst.copy()
        return fstlib.core.Fst(newfst)

    def decode(self, mapper):
        self.fst.decode(mapper)
        return self

    def delete_arcs(self, state, n=0):
        self.fst.delete_arcs(state, n)
        return self

    def delete_states(self, states=None):
        self.fst.delete_states(states)
        return self

    def draw(self, filename, isymbols=None, osymbols=None, ssymbols=None, acceptor=False, title="", 
          width=8.5, height=11, portrait=False, vertical=False, ranksep=0.4, nodesep=0.25, 
          fontsize=14, precision=5, float_format='g', show_weight_one=False):
        
        self.fst.draw(filename, isymbols, osymbols, ssymbols, acceptor, title, width, height, portrait, 
                vertical, ranksep, nodesep, fontsize, precision, float_format, show_weight_one)

    def encode(self, mapper):
        self.fst.encode(mapper)
        return self

    def final(self, state):
        return self.fst.final(state)

    def fst_type(self):
        return self.fst.fst_type()

    def info(self, name=None):
        df = fstlib.ext.info(self, name=name)
        return df

    def input_symbols(self):
        return self.fst.input_symbols()

    def invert(self):
        self.fst.invert()
        return self

    def is_final(self, state):
        return self.final(state) != pywrapfst.Weight.zero(self.weight_type())

    def minimize(self, delta=fstlib.SHORTEST_DELTA, allow_nondet=False):
        self.fst.minimize(delta, allow_nondet)
        return self

    def mutable_arcs(self, state):
        return self.fst.mutable_arcs(state)

    def mutable_input_symbols(self):
        return self.fst.mutable_input_symbols()

    def mutable_output_symbols(self):
        return self.fst.mutable_output_symbols()

    def num_arcs(self, state):
        return self.fst.num_arcs(state)

    def num_input_epsilons(self, state):
        return self.fst.num_input_epsilons(state)

    def num_output_epsilons(self, state):
        return self.fst.num_output_epsilons(state)
    
    def num_paths(self):
        """ Returns the number of paths of the fst via shortest distance. 
        An extension to pywrapfst.Fst. """
        return fstlib.paths.get_number_of_paths_from_fst(self)

    def num_states(self):
        return self.fst.num_states()
    
    def num_states_if_known(self):
        return self.fst.num_states_if_known()

    def output_symbols(self):
        return self.fst.output_symbols()

    def paths(self):
        """ Runs a path depth first search and returns a list of paths.
        An extension to pywrapfst.Fst. """
        return fstlib.paths.get_paths_from_fst(self)

    def print(self, isymbols=None, osymbols=None, ssymbols=None, acceptor=False, show_weight_one=False, missing_sym=""):
        return self.fst.print(isymbols, osymbols, ssymbols, acceptor, show_weight_one, missing_sym)

    def project(self, project_type):
        self.fst.project(project_type)
        return self

    def properties(self, mask, test):
        return self.fst.properties(mask,test)

    def prune(self, delta=fstlib.DELTA, nstate=fstlib.NO_STATE_ID, weight=None):
        self.fst.prune(delta, nstate, weight)
        return self

    def push(self, delta=fstlib.SHORTEST_DELTA, remove_total_weight=False, reweight_type="to_initial"):
        self.fst.push(delta, remove_total_weight, reweight_type)
        return self

    def relabel_pairs(self, ipairs=None, opairs=None):
        self.fst.relabel_pairs(ipairs, opairs)
        return self

    def relabel_tables(self, old_isymbols=None, new_isymbols=None, unknown_isymbols="", attach_new_isymbols=True,
                    old_osymbols=None, new_osymbols=None, unknown_osymbols="", attach_new_osymbols=True):
        self.fst.relabel_tables(old_isymbols, new_isymbols, unknown_isymbols, attach_new_isymbols,
                          old_osymbols, new_osymbols, unknown_osymbols, attach_new_osymbols)
        return self

    def reserve_arcs(self, state, n):
        return self.fst.reserve_arcs(state, n)

    def reserve_states(self, n):
        return self.fst.reserve_states(n)

    def reweight(self, potentials, reweight_type="to_initial"):
        self.fst.reweight(potentials, reweight_type)
        return self

    def rmepsilon(self, queue_type='auto', connect=True, weight=None, nstate=fstlib.NO_STATE_ID, delta=fstlib.SHORTEST_DELTA):
        self.fst.rmepsilon(queue_type, connect, weight, nstate, delta)
        return self

    def set_final(self, state, weight=None):
        self.fst.set_final(state, weight)
        return self

    def set_input_symbols(self, syms):
        self.fst.set_input_symbols(syms)
        return self

    def set_output_symbols(self, syms):
        self.fst.set_output_symbols(syms)
        return self

    def set_properties(self, props, mask):
        return self.fst.set_properties(props, mask)

    def set_start(self, state):
        self.fst.set_start(state)
        return self

    def start(self):
        return self.fst.start()

    def states(self):
        return self.fst.states()

    def to_dataframe(self, select='arcs', to_real=False):
        """ Converts an FST to a pandas DataFrame.
        An extension to pywrapfst.Fst. """
        if self.input_symbols() is not None:
            isyms = dict([(i, s) for i,s in self.input_symbols()])
        else:
            isyms = None
        if self.output_symbols() is not None:
            osyms = dict([(i, s) for i,s in self.output_symbols()])
        else:
            osyms = None

        def get_arc_record(state, arc):
            if isyms is not None:
                isym = isyms[arc.ilabel]
            else:
                isym = arc.ilabel
            if osyms is not None:
                osym = osyms[arc.olabel]
            else:
                osym = arc.olabel
            if to_real:
                weight = np.exp(-float(arc.weight))
            else:
                weight = float(arc.weight)
            record = (str(id(arc)), state, arc.nextstate, isym, osym, weight)
            return record

        if select == 'arcs':
            tuples = [get_arc_record(state, arc) for state in self.states() for arc in self.arcs(state)]
            df = pd.DataFrame.from_records(tuples, columns=['arcid', 'state_from', 'state_to', 'ilabel', 'olabel', 'weight'])
            df.set_index('arcid', inplace=True)
        elif select == 'states':
            tuples = [(state, np.exp(-float(self.final(state))) if to_real else float(self.final(state)), self.is_final(state)) for state in self.states()]
            index = [str(id(state)) for state in self.states()]
            df = pd.DataFrame.from_records(tuples, columns=['state', 'weight', 'is_final'], index=index)
        else:
            raise FSTError("Unknown select parameter, can be 'arcs' or 'states'")
        return df

    def to_graphviz(self, isymbols=None, osymbols=None, ssymbols=None, acceptor=False, title="", 
          width=8.5, height=11, portrait=True, vertical=False, ranksep=0.4, nodesep=0.25, 
          fontsize=14, precision=5, float_format='g', show_weight_one=False):
        import graphviz as gv
        
        tmpfile = tempfile.mktemp()
        self.draw(tmpfile, isymbols, osymbols, ssymbols, acceptor, title, width, height, portrait, 
                vertical, ranksep, nodesep, fontsize, precision, float_format, show_weight_one)
        with open(tmpfile, 'r') as fd:
            dot = fd.read()
        os.unlink(tmpfile)
        src = gv.Source(dot)
        return src

    def to_real(self):
        """ Shortcut for a constructive weight map from log to real, mainly
        for plotting and printing purposes. An extension to pywrapfst.Fst."""
        return fstlib.weight_map(self, fstlib.tools.neglog_to_real)

    def to_svg(self, **kwargs):
        g = self.to_graphviz(**kwargs)
        g.format='svg'
        return g

    def to_png(self, **kwargs):
        g = self.to_graphviz(**kwargs)
        g.format = 'png'
        return g

    def topsort(self):
        self.fst.topsort()
        return self

    def union(self, ifst):
        self.fst.union(ifst.fst)
        return self

    def verify(self):
        return self.fst.verify()

    def weight_map(self, func, with_final_states=True):
        """ Maps each arc and final state weight to a new value using the function provided.
         This is an extension to the pywrapfst Fst class. """
        for state in self.states():
            final_weight = self.final(state)
            if final_weight != pywrapfst.Weight.zero(self.weight_type()) and with_final_states: ## final state
                self.set_final(state, func(final_weight))
            mai = self.mutable_arcs(state)
            for arc in mai:
                arc.weight = pywrapfst.Weight(self.weight_type(), func(arc.weight))
                mai.set_value(arc)

        return self

    def weight_type(self):
        return self.fst.weight_type()

    def write(self, filename):
        return self.fst.write(filename)

    def write_to_string(self):
        return self.fst.write_to_string()

    @classmethod
    def read(cls, source):
        ofst = pywrapfst.Fst.read(source)
        return Fst(ofst)

    @classmethod
    def read_from_string(cls, state):
        ofst = pywrapfst.Fst.read_from_string(state)
        return Fst(ofst)

    def __getstate__(self):
        return self.write_to_string()

    def __setstate__(self, state):
        self.fst = pywrapfst.Fst.read_from_string(state)

class FSTError(Exception):
    def __init__(self, message="<unknown error>"):
        self.message = message

    def __str__(self):
        return str(self.message)

## function definitions
_last_seed = None

## function definitions
def generate_seed():
    global _last_seed
    random_data = os.urandom(4) 
    seed = int.from_bytes(random_data, byteorder="big") 
    _last_seed = seed
    return seed

def arcmap(ifst, delta=fstlib.DELTA, map_type='identity', power=1.0, weight=None):
    newfst = pywrapfst.arcmap(ifst.fst, delta, map_type, power, weight)
    return Fst(newfst)

def compact_symbol_table(syms):
    return pywrapfst.compact_symbol_table(syms)

def compose(ifst1, ifst2, compose_filter='auto', connect=True):
    newfst = pywrapfst.compose(ifst1.fst, ifst2.fst, compose_filter, connect)
    return Fst(newfst)

def convert(ifst, fst_type=""):
    newfst = pywrapfst.convert(ifst.fst, fst_type)
    return Fst(newfst)

def determinize(ifst, delta=fstlib.SHORTEST_DELTA, det_type='functional', nstate=fstlib.NO_STATE_ID, 
                subsequential_label=0, weight=None, incremental_subsequential_label=False):
    newfst = pywrapfst.determinize(ifst.fst, delta, det_type, nstate, subsequential_label, weight, incremental_subsequential_label)
    return Fst(newfst)

def difference(ifst1, ifst2, compose_filter='auto', connect=True):
    newfst = pywrapfst.difference(ifst1.fst, ifst2.fst, compose_filter, connect)
    return Fst(newfst)

def disambiguate(ifst, delta=fstlib.DELTA, nstate=fstlib.NO_STATE_ID, subsequential_label=0, weight=None):
    newfst = pywrapfst.disambiguate(ifst.fst, delta, nstate, subsequential_label, weight)
    return Fst(newfst)

# def divide(lhs, rhs):
#     return pywrapfst.divide(lhs, rhs)
# commented out to avoid forwarding overhead

def epsnormalize(ifst, eps_norm_type="input"):
    newfst = pywrapfst.epsnormalize(ifst.fst, eps_norm_type)
    return Fst(newfst)

def equal(ifst1, ifst2, delta=fstlib.DELTA):
    return pywrapfst.equal(ifst1.fst, ifst2.fst, delta)

def equivalent(ifst1, ifst2, delta=fstlib.DELTA):
    return pywrapfst.equivalent(ifst1.fst, ifst2.fst, delta)

def intersect(ifst1, ifst2, compose_filter='auto', connect=True):
    newfst = pywrapfst.intersect(ifst1.fst, ifst2.fst, compose_filter, connect)
    return Fst(newfst)

def isomorphic(ifst1, ifst2, delta=fstlib.DELTA):
    return pywrapfst.isomorphic(ifst1.fst, ifst2.fst, delta)

# def plus(lhs, rhs):
#     return pywrapfst.plus(lhs, rhs)
# commented out to avoid forwarding overhead

# def power(lhs, rhs):
#     return pywrapfst.power(lhs, rhs)
# commented out to avoid forwarding overhead

def prune(ifst, delta=fstlib.DELTA, nstate=fstlib.NO_STATE_ID, weight=None):
    newfst = pywrapfst.prune(ifst.fst, delta, nstate, weight)
    return Fst(newfst)

def push(ifst, delta=fstlib.DELTA, push_weights=False, push_labels=False, 
         remove_common_affix=False, remove_total_weight=False, reweight_type="to_initial"):
    newfst = pywrapfst.push(ifst.fst, delta, push_weights, push_labels, remove_common_affix, remove_total_weight, reweight_type)
    newfst = Fst(newfst)
    if remove_total_weight:
        logger.warn("'remove_total_weight' is not working in pywrapfst. Use destructive method instead. We're using a workaround here which might or might not do what you expect.")
        # for state in newfst.states():
        #     if newfst.is_final(state):
        #         newfst.set_final(state, fstlib.Weight.one(newfst.weight_type()))
    return newfst

def randequivalent(ifst1, ifst2, npath=1, delta=fstlib.DELTA, select='uniform', max_length=fstlib.MAX_INT32, seed=None):
    if seed is None:
        seed = generate_seed()
    return pywrapfst.randequivalent(ifst1.fst, ifst2.fst, npath, delta, select, max_length, seed)

def randgen(ifst, npath=1, select='uniform', max_length=fstlib.MAX_INT32, weighted=False, remove_total_weight=False, seed=None):
    if seed is None:
        seed = generate_seed()
    newfst = pywrapfst.randgen(ifst.fst, npath, select, max_length, weighted, remove_total_weight, seed)
    return Fst(newfst)

def replace(pairs, call_arc_labeling='input',  return_arc_labeling='neither', 
            epsilon_on_replace=False, return_label=0):
    return pywrapfst.replace(pairs, call_arc_labeling, return_arc_labeling, epsilon_on_replace, return_label)

def reverse(ifst, require_superinitial=True):
    newfst = pywrapfst.reverse(ifst.fst, require_superinitial)
    return Fst(newfst)

def shortestdistance(fst, delta=fstlib.SHORTEST_DELTA, nstate=fstlib.NO_STATE_ID, queue_type='auto', reverse=False):
    return pywrapfst.shortestdistance(fst.fst, delta, nstate, queue_type, reverse)

def shortestpath(fst, delta=fstlib.SHORTEST_DELTA, nshortest=1, nstate=fstlib.NO_STATE_ID, queue_type='auto', unique=False, weight=None):
    newfst = pywrapfst.shortestpath(fst.fst, delta, nshortest, nstate, queue_type, unique, weight)
    return Fst(newfst)

def statemap(ifst, map_type):
    newfst = pywrapfst.statemap(ifst.fst, map_type)
    return Fst(newfst)

def synchronize(ifst):
    newfst = pywrapfst.synchronize(ifst.fst)
    return Fst(newfst)

# def times(lhs, rhs):
#     return pywrapfst.times(lhs, rhs)
# commented out to avoid forwarding overhead
