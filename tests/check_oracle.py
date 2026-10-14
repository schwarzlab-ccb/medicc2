import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
"""Independent min-plus oracle, incremental/full checks, and real-seed pilots."""
import copy
import csv
import gc
import itertools
import json
import math
import sys
import time
from io import StringIO
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from Bio import Phylo
from medicc.full_fst import ancestors as hp


from oracle import get_parts, Primitive, relation, minplus, optimum


def main():
    started = time.perf_counter()
    rng = np.random.default_rng(20261013)
    parts, symbols = get_parts()
    counts = dict(oracle_messages=0, oracle_profiles=0, directed_edges=0,
                  incremental_neighbors=0, shared_cache_checks=0, parallel_scores=0,
                  nni_searches=0, metric_checks=0)
    for template in ['##', '#X#']:
        words = [''.join(next(it) if c=='#' else c for c in template)
                 for digits in itertools.product('012345678',repeat=2) for it in [iter(digits)]]
        ix = {s:i for i,s in enumerate(words)}
        primitive = {k:relation(Primitive(parts[k]),words,symbols) for k in ['LOH','LG','W1','L_LOH','G']}
        # Literal full two-WGD transducer relation on ALL intermediate strings.
        a = minplus(primitive['L_LOH'],primitive['G'])
        t = a
        for _ in range(2):
            t = minplus(minplus(t,primitive['W1']),a)

        def oracle(tree, samples, profiles, cache):
            expected = {k:{v:0} for k,v in samples.items()}
            for node in tree.find_clades(order='postorder'):
                if node.clades and node is not tree.root:
                    a,b = node.clades
                    expected[node.name] = optimum(words,t,expected[a.name],expected[b.name])
                    assert cache[node.name].native.language() == expected[node.name]
                    counts['oracle_messages'] += 1
            for node in tree.find_clades():
                for child in node.clades:
                    p = profiles[node.name]; q = profiles[child.name]
                    exact = t[ix[p],ix[q]]
                    assert hp.model().scorer.score(p,q,directed=True) == exact
                    assert hp.model().scorer.score(*hp.Model.shorten(p,q),directed=True) == exact
                    counts['directed_edges'] += 1
                    if child.clades:
                        candidates = expected[child.name]
                        value = min(t[ix[p],ix[w]]+c for w,c in candidates.items())
                        assert t[ix[p],ix[q]]+candidates[q] == value
                        counts['oracle_profiles'] += 1

        for case in range(8):
            tree = Phylo.read(StringIO('(((a,b)v,(c,d)w)r)internal_0;'),'newick')
            samples = {k:hp.Genome(words[int(rng.integers(81))]) for k in 'abcd'}
            samples['internal_0'] = hp.Genome('1X1' if 'X' in template else '11')
            prof, cache = hp.reconstruct(tree,samples,hp.FULL_TWO,hp.FULL_TWO,'internal_0',upper_cache=True,n_cores=8)
            oracle(tree,samples,prof,cache)
            before = {k:v.native.language() for k,v in cache.items()}
            neighbors = list(hp.medicc.nni.nni_neighbors(tree))
            reference = {}
            for neighbor, move in neighbors:
                inc, ic = hp.reconstruct(neighbor,samples,hp.FULL_TWO,hp.FULL_TWO,'internal_0',old_uppass_cache=cache,nni_move=move)
                full, fc = hp.reconstruct(neighbor,samples,hp.FULL_TWO,hp.FULL_TWO,'internal_0',upper_cache=True)
                assert inc == full
                assert {k:v.native.language() for k,v in ic.items()} == {k:v.native.language() for k,v in fc.items()}
                oracle(neighbor,samples,inc,ic)
                hp.update_lengths(neighbor,hp.FULL_TWO,inc)
                reference[hp.medicc.tree_hash.get_topology_hash(neighbor)] = hp.medicc.tools.sum_of_branch_length(neighbor)
                counts['incremental_neighbors'] += 1
            results = hp.medicc.nni.evaluate_nni_neighbors_parallel(tree,cache,samples,'internal_0',0,8,hp.FULL_TWO,hp.FULL_TWO,{})
            for neighbor, prof2, _, score, _, _ in results:
                assert reference[hp.medicc.tree_hash.get_topology_hash(neighbor)] == score
                counts['parallel_scores'] += 1
            assert before == {k:v.native.language() for k,v in cache.items()}
            counts['shared_cache_checks'] += 1
            if case == 0:
                a = hp.medicc.core.nni_mode(copy.deepcopy(tree),samples,hp.FULL_TWO,hp.FULL_TWO,normal_name='internal_0',n_cores=1,nni_stop_on_plateau=True,nni_max_frontier=1,nni_max_iter=math.inf)
                b = hp.medicc.core.nni_mode(copy.deepcopy(tree),samples,hp.FULL_TWO,hp.FULL_TWO,normal_name='internal_0',n_cores=8,nni_stop_on_plateau=True,nni_max_frontier=1,nni_max_iter=math.inf)
                assert a['best_score']==b['best_score'] and a['trace']==b['trace']
                assert hp.medicc.tree_hash.get_topology_hash(a['best_trees'][0]) == hp.medicc.tree_hash.get_topology_hash(b['best_trees'][0])
                counts['nni_searches'] += 2
            del cache, ic, fc
            gc.collect()
        print(template,counts,flush=True)

    # Distinguish the two-WGD library from the one/three-WGD pair libraries.
    assert hp.model().scorer.score('1X1','8X8',directed=True)==4
    summary=dict(status='ok', mismatches=0, counts=counts,seconds=time.perf_counter()-started)
    print(json.dumps(summary,indent=2),flush=True)
    return summary

if __name__=='__main__':
    main()
