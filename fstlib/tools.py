import numpy as np
import pandas as pd
import functools
import logging
import fstlib
import fstlib.core
import fstlib.algos
import fstlib.paths

logger = logging.getLogger(__name__)

class FstToolsError(Exception):
    pass

def strings(infst, delimiter="", to_real=False):
    """ returns all strings from a fsa """
    ipaths=list()
    opaths=list()
    scores=list()

    isyms = dict([(i,s if isinstance(s, str) else str(s, 'utf-8')) for i,s in infst.input_symbols()])
    osyms = dict([(i,s if isinstance(s, str) else str(s, 'utf-8')) for i,s in infst.output_symbols()])
    for path in fstlib.paths.get_paths_from_fst(infst):
        iseq = delimiter.join([isyms[p.ilabel] if isyms is not None else str(p.ilabel) for p in path])
        oseq = delimiter.join([osyms[p.olabel] if osyms is not None else str(p.olabel) for p in path])
        if infst.arc_type() == fstlib.Semiring.LOG or infst.arc_type() == fstlib.Semiring.TROPICAL:
            score = functools.reduce(fstlib.times, [p.weight for p in path] + [path.finalWeight,])
            #score = np.sum([float(p.weight) for p in path] + [float(path.finalWeight),])
        else:
            raise FstToolsError('semiring not implemented')
        ipaths.append(iseq)
        opaths.append(oseq)
        scores.append(float(score))

    digits = int(np.log10(len(ipaths)))+1
    result = pd.DataFrame(data=list(zip(ipaths, opaths, scores)), columns=['input', 'output','weight'])
    result.sort_values('weight', inplace=True)
    result.index=[("path%%.%dd" % digits) % i for i in range(len(ipaths))]
    if to_real:
        result.weight = np.exp(-result.weight)

    return result

def mass_intersect_quick(ifsts, prune_nstate = -1, prune_weight = "", prune_level = 0, determinize=True, minimize=True, rmepsilon=True, delta_det=fstlib.SHORTEST_DELTA, delta_min=fstlib.SHORTEST_DELTA, return_filename=False):
    return multicommand(fstlib.core.intersect, ifsts, prune_nstate, prune_weight, prune_level, determinize=determinize, 
                            minimize=minimize,
                          rmepsilon=rmepsilon, sort='olabel', return_filename=return_filename, delta_min=delta_min, 
                          delta_det=delta_det)

def multicommand(func, ifsts, prune_nstate = -1, prune_weight = "", prune_level = 0, determinize=True, minimize=True, rmepsilon=True,
                    sort=None, repeat=1, delta_det=fstlib.SHORTEST_DELTA, delta_min=fstlib.SHORTEST_DELTA, return_filename=False):
    """ applies a function multiple times to multiple fsts or repeatedly to the same.
    Uses Cyril's "log" strategy to speed up computation. """
        
    logger.info("Multicommand %s" % str(func))

    remaining = ifsts[:]
    level = 0
    funcalls = 0
    ofst = None
    while len(remaining)>1:
        next_level = remaining[:]
        for i in range(len(remaining)//2): ## build pairs
            ind1 = 2*i
            ind2 = 2*i+1
            ifst1 = remaining[ind1]
            ifst2 = remaining[ind2]
                
            next_level.remove(ifst1)
            next_level.remove(ifst2)

            #logger.info("processing %d and %d." % (ind1, ind2))	
            if sort == 'olabel': ## sort ifst1
                ifst1.arcsort(sort_type=sort)
            elif sort == 'ilabel': ## sort ifst2
                ifst2.arcsort(sort_type=sort)
            else: ## do not sort
                pass
            
            ## do the function
            ofst = func(ifst1, ifst2)
            funcalls += 1

            do_prune = False
            if (prune_weight != "" or prune_nstate >=0) and prune_level <= level:
                do_prune = True

            if not determinize and do_prune: ## if we are determinizing use weighted determinization instead
                ofst = fstlib.core.prune(ofst, weight=prune_weight, nstate=prune_nstate)

            if determinize:
                if rmepsilon:
                    ofst.rmepsilon()
                if do_prune:
                    ofst = fstlib.core.determinize(ofst, delta=delta_det, nstate=prune_nstate, weight=prune_weight)
                else:
                    ofst = fstlib.core.determinize(ofst)
            if minimize:
                ofst.minimize(delta = delta_min)							
            
            next_level.append(ofst)
        remaining = next_level
        level += 1
        
    logger.info("done (%d function calls)!" % funcalls)
        
    return ofst

def _exec_func(func, ifst1, ifst2, delta_min, output, determinize=False, minimize=False, rmepsilon=False, sort='olabel'):
    if sort == 'olabel': ## sort ifst1
        ifst1.arcsort(sort_type=sort)
    elif sort == 'ilabel': ## sort ifst2
        ifst2.arcsort(sort_type=sort)
    else: ## do not sort
        pass
            
    ## do the function
    ofst = func(ifst1, ifst2)

    if determinize:
        if rmepsilon:
            ofst.rmepsilon()
        ofst = fstlib.core.determinize(ofst)
    if minimize:
        ofst.minimize(delta = delta_min)
    output.put(ofst)


def strings_to_count_matrix(list_of_strings, symbols=None):
    """ transforms a list of strings into a count matrix by counting the number of occurences of 
    each symbol in all strings at each position in the list. """
    if symbols is None:
        symbols = list(set(''.join(list_of_strings)))
        symbols.sort()
    n = np.max([len(s) for s in list_of_strings]) ## determine length of longest strings
    m = len(symbols)
    count_matrix = np.zeros((m, n))
    for i in range(n):
        for j in range(len(list_of_strings)):
            s = list_of_strings[j][i]
            k = symbols.index(s)
            count_matrix[k,i] += 1

    return count_matrix, symbols


def neglog_to_real(x):
    return np.exp(-float(x))

def real_to_neglog(x):
    return -np.log(float(x))
