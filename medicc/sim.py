import logging
import random

import Bio
import numpy as np
import pandas as pd
import scipy as sp
import copy

from medicc import tools

logger = logging.getLogger(__name__)

def evolve(cnstr, mu=5, plen=0.2, pwgd=0.05, pgain=0.5, maxcn=8, mincn=0, verbose=True, seed=None, maxwgd=3):
    """ Evolves a copy-number string, e.g. 11111X11111 with the given parameters. """
    nevents = sp.stats.poisson.rvs(mu, size=1)[0]
    cnstr_new = np.array(list(cnstr))
    nloss = ngain = nwgd = 0
    if seed is not None:
        np.random.seed(seed)
    for i in range(nevents):
        if nwgd < maxwgd and np.random.uniform() <= pwgd: ## WGD
            nwgd += 1
            _mutate(cnstr_new, 'gain', 0, len(cnstr_new), maxcn, mincn, ignore_chr_boundaries=True, verbose=verbose)
        else:
            available_positions = np.flatnonzero((cnstr_new != '0') & (cnstr_new != 'X'))
            if len(available_positions) == 0:
                logger.warning("No available positions to mutate.")
                break
            start = np.random.choice(available_positions)
            length = np.random.geometric(plen)
            end = min(start+length, np.max(available_positions)+1)
            if np.random.uniform() <= pgain: ## gain
                ngain += 1
                _mutate(cnstr_new, 'gain', start, end, maxcn, mincn, verbose=verbose)
            else:
                nloss += 1
                _mutate(cnstr_new, 'loss', start, end, maxcn, mincn, verbose=verbose)
    return ("".join(cnstr_new), nloss, ngain, nwgd)

def _mutate(cnstr, event, start, end, maxcn, mincn, ignore_chr_boundaries=False, verbose=True):
    """ Mutates an individual copy-number string, e.g. 1111X111111. """
    if verbose:
        logger.info(f"{event} event from {start} to {end}")
    for i in range(start, end):
        if cnstr[i]=='X':
            if ignore_chr_boundaries:
                continue
            else:
                break
        oldcn = tools.hex2int(cnstr[i])
        if event == 'gain':
            if oldcn != 0:
                newcn = oldcn + 1
            else:
                newcn = oldcn
            if newcn > maxcn:
                newcn = maxcn
        elif event == 'loss':
            newcn = oldcn - 1
            if newcn < mincn:
                newcn = mincn
        else:
            logger.error('unknown event')
            newcn = oldcn
        cnstr[i] = tools.int2hex(newcn)

def evolve_along_tree(tree, normal_sequence, normal_name='diploid', rate_multiplier=1, 
    plen=0.2, pwgd=0.05, pgain=0.5, maxcn=8, mincn=0, maxwgd=3, 
    verbose=True, seed=None, inplace=True):
    """ Mutates one allele along a given tree starting from the normal sample with name indicated
    and normal sequence as given. All event rates are globally multiplied by rate_multiplier. """
    if not inplace:
        tree = copy.deepcopy(tree)
    stack = []
    child_counter = {}
    finished = False
    current_clade = tree.root
    child_counter[current_clade]=0
    for clade in tree.root.clades:
        if clade.name == normal_name: 
            clade.seq = normal_sequence
            diploid_bl = clade.branch_length
            current_clade.seq = evolve(normal_sequence, mu=diploid_bl*rate_multiplier, 
                plen=plen, pwgd=pwgd, pgain=pgain, maxcn=maxcn, mincn=mincn, maxwgd=maxwgd,
                verbose=verbose, seed=seed)[0]
            break

    while not finished:
        i = child_counter[current_clade]
        if i < len(current_clade.clades):
            if current_clade.clades[i].name == normal_name:
                child_counter[current_clade]+=1
            else:
                child_clade = current_clade.clades[i]
                child_clade.seq = evolve(current_clade.seq, mu=child_clade.branch_length*rate_multiplier, 
                    plen=plen, pwgd=pwgd, pgain=pgain, maxcn=maxcn, mincn=mincn, maxwgd=maxwgd,
                    verbose=verbose, seed=seed)[0]
                stack.append(current_clade)
                child_counter[current_clade]+=1
                current_clade = child_clade
                child_counter[current_clade]=0
        else:
            try:
                current_clade = stack.pop()
            except IndexError:
                finished = True

    return tree
    

def rcoal(n, tips = None, uniform_branch_length=False, normal_name='diploid'):
    """As in rcoal function in R package Ape, function generates a random tree and random branch lengths."""
    x = np.divide(2*np.random.exponential(size = n-1),np.multiply(np.array(list(reversed(range(2,n+1)))), np.array(list(reversed(range(1, n))))))
    if (n == 2):
        edge = pd.DataFrame()
        edge['a'] = [2,2]
        edge['b'] = [0,1]
        edge['name_a'] = pd.Series(['mrca','mrca'])
        if tips and len(tips)==n:
            edge['name_b'] = pd.Series([tips[0],tips[1]])
        else:
            edge['name_b'] = pd.Series(['sp_0001','sp_0002'])
        edge['edge_length'] = np.array([x,x])
        reordered = edge.copy()
    elif (n == 3):
        edge = pd.DataFrame()
        edge['a'] = [3,4,4,3]
        edge['b'] = [4,0,1,2]
        edge['name_a'] = pd.Series(['mrca','internal_1','internal_1','mrca'])
        if tips and len(tips)==n:
            edge['name_b'] = pd.Series(['internal_1',tips[0],tips[1],tips[2]])
        else:
            edge['name_b'] = pd.Series(['internal_1','sp_0001','sp_0002','sp_0003'])
        edge['edge_length'] = np.array([x[1],x[0],x[0],sum(x)])
        reordered = edge.copy()
    else:
        edge = np.zeros((2*n - 2,2), dtype=float)
        edge_length = np.array([0]*(2*n - 2),dtype = float)
        h = np.array([0]*(2*n - 1),dtype = float)
        node_height = np.cumsum(x)
        pool = np.array(list(range(0,n)))# 1:n
        nextnode  = 2*n - 2
        for i in range(0,n-1):
            y = random.sample(list(pool), 2)
            ind = np.array([i*2 , i*2 + 1])
            edge[ind, 1] = y
            edge[ind, 0] = nextnode
            edge_length[ind] = node_height[i] - h[y]
            h[nextnode] = node_height[i]
            pool = np.append(pool[np.isin(pool,y,invert = True)],[nextnode],axis = 0)
            nextnode -= 1
        edge = pd.DataFrame(edge,dtype = object)
        edge.columns = ['a','b']


        reordered = pd.DataFrame(columns=['a','b'])
        root = edge.index[-2]
        parent = [edge.iloc[root][0]] #stack
        child = edge.iloc[root][1]
        reordered = pd.concat([reordered, edge.iloc[[root]]])
        edge.drop(edge.index[root], inplace = True)
        while  not edge.empty:
            if (set (edge.a.isin([parent[-1]]))) == {False}:#checks if both occurences took place, and value can be popped out of stack
                del parent[-1]
            if child <= n :
                root = edge.a.isin([parent[-1]]).idxmax()
                child = int(edge.iloc[edge.index==root]['b'])
                reordered = pd.concat([reordered, edge.iloc[edge.index==root]])
                edge.drop(root, inplace = True)
                if child < n:
                    if (set (edge.a.isin([parent[-1]]))) == {False}:
                        del parent[-1]
            else:
                parent.append(child)
                root =  edge.a.isin([child]).idxmax()
                child = int(edge.iloc[edge.index==root]['b'])
                reordered = pd.concat([reordered, edge.iloc[edge.index==root]])
                edge.drop(root, inplace = True)
                if child < n:
                    if (set (edge.a.isin([parent[-1]]))) == {False}:
                        del parent[-1]
        reordered['a'] = reordered['a'].astype(int)
        reordered['b'] = reordered['b'].astype(int)
        reordered.loc[:,'edge_length'] = pd.Series(edge_length).reindex(index=reordered.index)
        reordered.replace((reordered['b'].loc[reordered['b']<n].values),list(range(0,n)), inplace = True)
        reordered.edge_length = round(reordered.edge_length,3)
        reordered.index = range(len(reordered))

        if tips and len(tips)==n:
            reordered['name_b'] = [tips[reordered.b[i]] if (reordered.b[i] < n)  else '' for i in reordered.index]
        else:
            reordered['name_b'] = ["sp_{0:04d}".format(reordered.b[i] +1) if (reordered.b[i] < n)  else '' for i in reordered.index]
        reordered.loc[reordered.name_b=='', 'name_b'] = ["internal_{0:d}".format(i +1) for i in range(len(reordered[reordered.name_b=='']))]
        internal_names =  reordered[reordered.b>n][['b','name_b']].copy()
        internal_names.set_index(['b'], inplace = True)
        internal_names = pd.concat([internal_names, pd.DataFrame(data = ['mrca'], index= [n], columns = ['name_b'])])
        reordered['name_a'] = [internal_names.name_b .loc[i]for i in reordered.a]

    reordered = reordered[['a','b','name_a', 'name_b', 'edge_length']]
    if uniform_branch_length:
        reordered['edge_length'] = 1
    return _edgelist_to_tree(reordered, add_normal=normal_name is not None, normal_name=normal_name)


def _edgelist_to_tree(el, add_normal=True, normal_name='diploid'):
    clades = {name:Bio.Phylo.PhyloXML.Clade(name=name) for name in el[['name_a','name_b']].stack(future_stack=True).unique()}
    if add_normal:
        clades[normal_name] = Bio.Phylo.PhyloXML.Clade(name=normal_name, branch_length=1)
        clades['mrca'].clades.append(clades[normal_name])

    for _, row in el.iterrows():
        parent = clades[row['name_a']]
        child = clades[row['name_b']]
        child.branch_length = row['edge_length']
        parent.clades.append(child)
    
    tree = Bio.Phylo.PhyloXML.Phylogeny(root = clades['mrca'])
    return tree
