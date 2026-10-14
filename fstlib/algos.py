'''
Created on 17 May 2010

@author: rfs

FST algorithms
'''

import array
import functools
import logging
import collections

import numpy as np

import fstlib.core
#from fstlib.core import *
import fstlib.tools

log = logging.getLogger(__name__)

class FstAlgoError(Exception):
    pass

class ConfusionNetworkFromLattice(object):
    """ Expects lattice in the form of a FST object. 
    Requires lattice to be topologically sorted
    (arcs are always directed from lower to higher state ids). It might require (not sure)
   the initial state to have index 0. """
    def __init__(self, lattice):
        self.lattice = lattice
        self.state2set = None ## state set ids for each state in lattice
        self.statesets = None ## list of lists of state ids contained in the state sets
        self.edgesets = None ## new edges between state sets
        
    def run(self, collapse_arcs=True):
        log.info('Constructing confusion network...')
        lattice = self.lattice
        outfst = fstlib.core.Fst(arc_type = self.lattice.arc_type())
        outfst.set_input_symbols(self.lattice.input_symbols())
        outfst.set_output_symbols(self.lattice.output_symbols())
        timeslices = fstlib.shortestdistance(fstlib.arcmap(fstlib.arcmap(fstlib.arcmap(self.lattice, map_type='rmweight'), map_type='times', weight=1), map_type='to_std'))
        statemap = np.array([int(float(t)) for t in timeslices])
        newstates = np.array(outfst.add_states(np.unique(statemap).shape[0]))
        outfst.set_start(newstates[0])
        for s in newstates:
            weights = [self.lattice.final(t) for t in np.flatnonzero(statemap==s) if self.lattice.is_final(t)]
            if len(weights)>0: ## newstate is final if any of the merged states were final (have non-inf weights)
                total_weight = functools.reduce(fstlib.times, weights)
                outfst.set_final(s, total_weight)

        for state in self.lattice.states():
            newstate = statemap[state]
            for arc in self.lattice.arcs(state):
                newnextstate = statemap[arc.nextstate]
                outfst.add_arc(newstate, (arc.ilabel, arc.olabel, arc.weight, newnextstate))

        if collapse_arcs:
            outfst=fstlib.determinize(outfst)
        return outfst

class PosteriorDecoding(object):
    """ expects pFSA in binary (i.e. openfst) format """
    def __init__(self, ifst):
        self.ifst = ifst

    def run(self):
        alpha = fstlib.shortestdistance(self.ifst)
        beta = fstlib.shortestdistance(self.ifst, reverse=True)
        
        outfst = fstlib.core.Fst(arc_type=self.ifst.arc_type())
        outfst.set_input_symbols(self.ifst.input_symbols())
        outfst.set_output_symbols(self.ifst.output_symbols())
        
        for state in self.ifst.states():
            newstate = outfst.add_state()
            if self.ifst.is_final(state):
                outfst.set_final(newstate, self.ifst.final(state))
            if state != newstate:
                raise FSTlibAlgosError('State mapping error') 
                ## this must not happen; I'm not sure it can currently due to the way the states are handled
        initial_state = self.ifst.start()
        outfst.set_start(initial_state)
        for state in self.ifst.states():
            for arc in self.ifst.arcs(state):
                posterior = float(alpha[state]) + float(arc.weight) + float(beta[arc.nextstate]) - float(beta[initial_state])
                outfst.add_arc(state, (arc.ilabel, arc.olabel, posterior, arc.nextstate))
        return outfst

class AbstractCountDP(object):
    def __init__(self, fst, trans_cache=None, state_depth_cache=None):
        """ constructor """
        self.fst = fst
        
        ## if transition cache was given use it, otherwise create it
        if trans_cache == None:
            self._build_trans_cache()
        else:
            self._trans_cache = trans_cache
        
        ## if state depth cache was given use it, otherwise create it
        if state_depth_cache == None:
            self._build_state_depth_cache()
        else:
            self._state_depth_cache = state_depth_cache
            
        assert(fst != None)

    def _build_trans_cache(self):
        #self._trans_cache = numpy.empty( (len(self.fst.states), len(self.fst.states)), dtype="object")
        self._trans_cache = [0] * len(self.fst.states)
        self._trans_cache = [[i] * len(self.fst.states) for i in self._trans_cache]
        
        for i in range(0, len(self.fst.states)):
            from_transs = [t for t in self.fst.arcs if t.state_from.id == i]
            for j in range(0, len(self.fst.states)):
                transs = len([t for t in from_transs if t.state_to.id == j])
                self._trans_cache[i][j] = transs
                
    def _build_state_depth_cache(self):
        """ tries to speed up the algo """
        self._state_depth_cache=list()
        depth = 0
        self._build_state_depth_cache_recurse(self.fst.states[0], depth)
        
    def _build_state_depth_cache_recurse(self, state, depth):
        """ recursion function of the pre inspection """
        
        if len(self._state_depth_cache) <= depth:
            self._state_depth_cache.append(set())
        self._state_depth_cache[depth].add(state.id)
        
        targets = list(set([t.state_to for t in state.get_outgoing_transitions()]))
        for s in targets:
            self._build_state_depth_cache_recurse(s, depth + 1)

    def max_seq_len(self):
        return len(self._state_depth_cache)-1

    def get_transition_cache(self):
        return self._trans_cache
    
    def get_state_depth_cache(self):
        return self._state_depth_cache
                
        
class ForwardCountDP(AbstractCountDP):
    """ Dynamic programming forward algorithm.
    This is not in log space and very generic. """
    
    def __init__(self, fst, trans_cache=None, state_depth_cache=None):
        """ constructor """
        super(ForwardCountDP, self).__init__(fst, trans_cache, state_depth_cache)
        self.alpha = None
        self._run()
        
    def _run(self):
        """ main algorithm """
        
        nstates = len(self.fst.states)
        ## initialize the forward matrix
        # the sequence index starts with 1, 0 means "no part of the sequence seen so far"
        ##self.alpha = numpy.zeros((nstates, 1))
        self.alpha = list()
        for i in range(0, self.max_seq_len() + 1):
            self.alpha.append(array.array("L", [0] * nstates))
        
        # we assume the starting state is always state 0
        self.alpha[0][0] = 1 

        for seqpos in range(1, self.max_seq_len() + 1):
            for stateid_to in self._state_depth_cache[seqpos]:
                value = 0
                for stateid_from in self._state_depth_cache[seqpos-1]:
                    value += self._trans_cache[stateid_from][stateid_to] * self.alpha[seqpos-1][stateid_from]
                self.alpha[seqpos][stateid_to] = value
    
class BackwardCountDP(AbstractCountDP):
    def __init__(self, fst, trans_cache=None, state_depth_cache=None):
        """ constructor """
        super(BackwardCountDP, self).__init__(fst, trans_cache, state_depth_cache)
        self.beta = None
        self._run()

    
    def _run(self):
        """ main algorithm """
        
        ## initialize the forward matrix
        # the sequence index starts with 1, 0 means "no part of the sequence seen so far"
        nstates = len(self.fst.states)
        self.beta = list()
        for stateid_from in range(0, self.max_seq_len() + 1):
            self.beta.append(array.array("L", [0] * nstates))
        
        # all final states are assigned 1
        for id in self._state_depth_cache[self.max_seq_len()]:
            self.beta[self.max_seq_len()][id] = 1 

        for seqpos in range(self.max_seq_len()-1, -1, -1):
            for stateid_from in self._state_depth_cache[seqpos]:
                value = 0
                for stateid_to in self._state_depth_cache[seqpos + 1]:
                    value += self.beta[seqpos+1][stateid_to] * self._trans_cache[stateid_from][stateid_to]  
                self.beta[seqpos][stateid_from] = value

class PosteriorDecodingCountDP(AbstractCountDP):
    """ Dynamic programming algorithm for posterior decoding of counts """
    
    def __init__(self, fst):
        
        """ constructor """
        self.forward_algo = ForwardCountDP(fst)
        self.backward_algo = BackwardCountDP(fst,  
                                   self.forward_algo.get_transition_cache(),
                                   self.forward_algo.get_state_depth_cache())
        super(PosteriorDecodingCountDP, self).__init__(fst, 
                                                  self.forward_algo.get_transition_cache(),
                                                  self.forward_algo.get_state_depth_cache())        
        self.posterior = [dict() for i in range(0, self.max_seq_len())]
        self._run()
        
    def _run(self):
        """ main algorithm """
        nstates = len(self.fst.states)
        
        for pos in range(0, self.max_seq_len()):
            for state_from in self._state_depth_cache[pos]:
                all_transitions = self.fst.states[state_from].get_outgoing_transitions()
                for state_to in self._state_depth_cache[pos+1]:
                    transitions = [t for t in all_transitions if t.state_to.id == state_to]
                    fw = self.forward_algo.alpha[pos][state_from]
                    bw = self.backward_algo.beta[pos+1][state_to]
                    for trans in transitions:
                        try:
                            value = self.posterior[pos][trans.isymbol]
                        except KeyError:
                            value = 0
                        value += fw * np.exp(-trans.weight) * bw
                        self.posterior[pos][trans.isymbol] = value



class AbstractDP(object):
    def __init__(self, fst, trans_cache=None):
        """ constructor """
        self.fst = fst
        if trans_cache == None:
            self._build_cache()
        else:
            self._trans_cache = trans_cache
        assert(fst != None)

    def _get_transitions(self, state_from, state_to):
        """ returns all transitions between states with caching.
        Expects numerical state ids. """
        
        return self._trans_cache[state_from][state_to]
    
    def _build_cache(self):
        #self._trans_cache = numpy.empty( (len(self.fst.states), len(self.fst.states)), dtype="object")
        self._trans_cache = [None] * len(self.fst.states)
        self._trans_cache = [[i] * len(self.fst.states) for i in self._trans_cache]
        
        for i in range(0, len(self.fst.states)):
            from_transs = [t for t in self.fst.arcs if t.state_from.id == i]
            for j in range(0, len(self.fst.states)):
                transs = array.array('f', [t.weight for t in from_transs if t.state_to.id == j])
                self._trans_cache[i][j] = transs
                
                
        
class PosteriorDecodingDP(AbstractDP):
    """ Dynamic programming algorithm for posterior decoding """
    
    def __init__(self, fst):
        """ constructor """
        super(PosteriorDecodingDP, self).__init__(fst)        
        self.forward_algo = ForwardDP(self.fst)
        self.backward_algo = BackwardDP(self.fst, 
                                   self.forward_algo.get_sequence_length(), 
                                   self.forward_algo.get_transition_cache())
        self._trans_cache = self.forward_algo.get_transition_cache()
        self.posterior = [dict() for i in range(0, self.forward_algo.get_sequence_length())]
        self._run()
        
    def _run(self):
        """ main algorithm """
        
        
        for pos in range(0, self.forward_algo.get_sequence_length()):
            for state_from in range(0, len(self.fst.states)):
                for state_to in range(0, len(self.fst.states)):
                    transitions = self._get_transitions(state_from, state_to)
                    fw = self.forward_algo.alpha[state_from, pos]
                    bw = self.backward_algo.beta[state_to, pos+1]
                    for trans in transitions:
                        try:
                            value = self.posterior[pos][trans.isymbol]
                        except KeyError:
                            value = 0
                        value += fw * np.exp(-trans.weight) * bw
                        self.posterior[pos][trans.isymbol] = value
                        
        
class ForwardDP(AbstractDP):
    """ Dynamic programming forward algorithm.
    This is not in log space and very generic. """
    
    def __init__(self, fst, trans_cache=None):
        """ constructor """
        super(ForwardDP, self).__init__(fst, trans_cache)
        self.alpha = None        
        self._run()
        
    def _run(self):
        """ main algorithm """
        
        nstates = len(self.fst.states)
        ## initialize the forward matrix
        # the sequence index starts with 1, 0 means "no part of the sequence seen so far"
        ##self.alpha = numpy.zeros((nstates, 1))
        self.alpha = list()
        self.alpha.append(array.array("f", [0] * nstates))
        
        # we assume the starting state is always state 0
        self.alpha[0][0] = 1 

        j=0
        while j<len(self.alpha):
            j+=1
            column = array.array("f", [0] * nstates)
            any_weight = False        
            for i in range(0, nstates):
                value = 0
                for k in range(0, nstates):
                    transitions = self._trans_cache[k][i]
                    for t in transitions:
                        value += self.alpha[j-1][k] * np.exp(-t)  
                column[i] = value
                if value > 0.0: any_weight = True
            if any_weight:
                self.alpha.append(column)
    
    def get_sequence_length(self):
        return len(self.alpha)
    
    def get_transition_cache(self):
        return self._trans_cache
                
class BackwardDP(AbstractDP):
    """ Dynamic programming backward algorithm.
    This is not in log space and very generic. """
    
    def __init__(self, fst, seqlen, trans_cache=None):
        """ constructor """
        super(BackwardDP, self).__init__(fst, trans_cache)
        self.seqlen = seqlen
        self.beta = None
        self._run()
        
    def _run(self):
        """ main algorithm """
        
        ## initialize the forward matrix
        # the sequence index starts with 1, 0 means "no part of the sequence seen so far"
        self.beta = np.zeros((len(self.fst.states), self.seqlen + 1))
        
        # all final states are assigned 1
        for state in self.fst.states:
            if state.is_final:
                self.beta[state.id, self.seqlen] = 1 

        for j in range(self.seqlen-1, -1, -1):
            for i in range(0, len(self.fst.states)):
                value = 0
                for k in range(0, len(self.fst.states)):
                    transitions = self._get_transitions(i, k)
                    for t in transitions:
                        value += self.beta[k, j+1] * np.exp(-t)  
                self.beta[i,j] = value
        
def posterior_decoding(ifst, collapse_arcs=True):
    pd = PosteriorDecoding(ifst).run()
    conf = ConfusionNetworkFromLattice(pd).run(collapse_arcs)
    return conf

class FSTlibAlgosError(Exception):
    pass
