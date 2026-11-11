import numpy as np
from scipy import stats


def _remove_duplicates(D):
    """ removes duplicates (entries with D=0) starting from the highest index. """
    finished=False
    while not finished:
        triu = np.triu_indices(D.shape[0], k=1)
        duplicates = np.flatnonzero(D[triu]==0)
        if len(duplicates)>0:
            dupe = duplicates[0]
            sel = max(triu[0][dupe], triu[1][dupe])
            D = np.delete(D, sel, axis=0)
            D = np.delete(D, sel, axis=1)
        else:
            finished=True
    return D

def star_topology_test(D):
    if D.shape[0] == 2:
        return 1
    """ Tests the distance matrix D for star topology. """
    ## check if any non-diagonal entry is zero
    D = _remove_duplicates(D)

    ## extract lower triangle matrix
    p=D[np.triu_indices(D.shape[0], k=1)]

    ## set up |p| x |samples| matrix T of zeros and ones mapping branch lengths
    ## in a star topology setting to the pairwise distances (e.g. dist(a,b) = branch_a + branch_b)
    count=0
    T = np.zeros((len(p), D.shape[1]))
    for i in range(T.shape[1]-1):
        for j in range(i+1, T.shape[1]):
            T[count,[i,j]]=1
            count = count + 1

    ## now find the optimal branch lengths b for our p
    Tinv = np.linalg.inv(np.matmul(T.T, T))
    bopt = np.matmul(np.matmul(p.T, T), Tinv).T
 
    ## scale the residuals by the square root of their variance (their sd)
    ## the estimates p are normal distributed with mean t (time of divergence) and variance t.
    res = (np.matmul(T, bopt)-p) / np.sqrt(p)
    res = np.matmul(res.T, res)

    ## the sum of squared residuals is therefore chisquare distributed with |p|-|samples| degrees
    ## of freedom. The "-|samples|" is because we optimize over the branch lengths (Tim Massingam,
    ## personal communication)
    
    prob = stats.chi2.cdf(res, len(p)-D.shape[1])
    retval = 1-prob
    return retval ## the p-value

def molecular_clock_test(D, normal_index=0):
    ## distances from diploid to all others but diploid
    p=np.delete(D, normal_index, axis=1)[normal_index, :] 
    ## remove other zeros
    p = p[p!=0]

    ## scale the residuals by the square root of their variance (their sd)
    ## the estimates p are normal distributed with mean t (time of divergence) and variance t.
    res = (np.mean(p) - p) / np.sqrt(p)
    res = np.matmul(res.T, res)
    
    ## the sum of squared residuals is therefore chisquare distributed with |p| degrees
    ## of freedom. (Tim Massingam, personal communication)

    prob = stats.chi2.cdf(res, len(p))
    return 1-prob ## the p-value


def MEDICCStatsError(Exception):
    pass
