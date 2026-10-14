import numpy as np
import pandas as pd
import logging
import fstlib
import fstlib.core
import fstlib.algos

logger = logging.getLogger(__name__)

class FstPathsError(Exception):
    pass

class Path(list):
    """ simple extension of list class that includes a final weight """
    __attributes = ["finalWeight"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args)
        self.finalWeight = kwargs.get("finalWeight", 0)
        valid = [p in self.__attributes for p in kwargs.keys()]
        if not all(valid):
            raise FstPathsError("Unknown attribute in [" + ",".join(kwargs.keys()) + "]. Allowed arguments are: [" + ",".join(self.__attributes) + "].")

    def copy(self):
        newpath = fstlib.Path(self)
        newpath.finalWeight = self.finalWeight
        return newpath

def get_paths_from_fst(ifst):
    """ generator that returns paths from a given fst.
    Doesn't work correctly if multiple final states along a single path, t.b.f. """

    if ifst.num_states()==0:
        return

    zero = fstlib.Weight.zero(ifst.weight_type())
    state = ifst.start()
    next_arc = {}
    stack = []
    path = Path()
    arc_cache = {}

    while state is not None:
        try:
            arcs = arc_cache[state]
        except KeyError:
            arcs = [a for a in ifst.arcs(state)]
            arc_cache[state] = arcs

        try:
            i = next_arc[state]
        except KeyError:
            i = 0
        
        #if self.fst.num_arcs(state) > i: ## move forward
        if len(arcs) > i: ## move forward
            arc = arcs[i]  ## get forward transition
            path.append(arc) ## save transition
            stack.append(state)
            next_arc[state] = i + 1 ## increase trans counter by 1
            state = arc.nextstate ## move to next state
        else: ## move back
            path.finalWeight = ifst.final(state)
            if path.finalWeight != zero:
                yield path.copy()

            next_arc[state] = 0
            try:
                state = stack.pop()
                path.pop()
            except IndexError:
                state = None

def get_number_of_paths_from_fst(ifst):
    """ Returns the total number of paths for a given fst. """
    if ifst.properties(fstlib.FstProperties.CYCLIC, True) == fstlib.FstProperties.CYCLIC:
        npaths = np.inf
    else:
        fsto = fstlib.arcmap(ifst, map_type='rmweight')
        fsto = fstlib.arcmap(fsto, map_type='to_log')
        sd = fstlib.shortestdistance(fsto, reverse=True)[0]
        npaths = int(np.round(np.exp(-float(sd))))
    return npaths
