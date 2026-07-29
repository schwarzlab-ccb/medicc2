import copy
import logging

import numpy as np
import pandas as pd
from Bio.Phylo.Consensus import _BitString, get_support

from medicc import io, tools
from medicc.core import main

# tqdm can be used for progress bars
try: 
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(x, disable):
        return x

logger = logging.getLogger(__name__)

def chr_wise_bootstrap_df(input_df):
    """Creates a bootstrap dataframe based on the original data. From the N original chromosomes in
    the dataset we draw N chromosomes *with* replacement to create a new bootstrap dataframe.
    """

    chroms = input_df.index.get_level_values('chrom').unique()
    cur_chroms = np.random.choice(chroms, len(chroms), replace=True)

    bootstrap_df = pd.DataFrame(columns=input_df.reset_index().columns)

    for i, cur_chrom in enumerate(cur_chroms):
        cur_data = input_df.reset_index().loc[(input_df.reset_index()['chrom'] == cur_chrom).values]
        cur_data['chrom'] = 'chr{}'.format(i+1)
        bootstrap_df = pd.concat([bootstrap_df, cur_data])

    bootstrap_df['chrom'] = tools.format_chromosomes(bootstrap_df['chrom'])
    bootstrap_df[['start', 'end']] = bootstrap_df[['start', 'end']].astype(int)
    bootstrap_df.set_index(['sample_id', 'chrom', 'start', 'end'], inplace=True)
    bootstrap_df.sort_index(inplace=True)

    logger.debug("Created chr-wise bootstrap dataframe")

    return bootstrap_df


def segment_wise_jacknife_df(input_df):
    """Creates a jackknife dataframe based on the original data. From the N original segments in
    the dataset we draw N segments *with* replacement and keeping their original order. For 
    duplicate segments only one copy will be kept. On average this will result in (1-1/e) * N segments. 
    """
    jacknife_df = input_df.copy()
    jacknife_df['segment'] = jacknife_df.reset_index()[['chrom', 'start']].astype(
        str).apply(lambda x: '-'.join(x), axis=1).values
    segments = jacknife_df['segment'].unique()

    cur_segments = np.unique(np.random.choice(segments, size=len(segments), replace=True))
    jacknife_df = jacknife_df.loc[jacknife_df['segment'].map(lambda x: x in cur_segments)].drop('segment', axis=1)

    logger.debug(
        "Created segment-wise jackknife dataframe with {} out of the {} segments".format(len(cur_segments), 
                                                                                         len(segments)))

    return jacknife_df


def bootstrap_shuffle_chroms(input_df, normal_name='diploid'):
    """Creates a bootstrap dataframe based on the original data. The chromosome-wise copy-number data 
    will be shuffled between samples. The corresponding dataframe will be generally similar to the
    original data but lack the distinct features.
    It can therefore not be used to calculate branch support or similar scores.
    """

    chroms = list(input_df.index.get_level_values('chrom').unique())
    sample_labels = [s for s in input_df.reset_index()['sample_id'].unique() if s != normal_name]
    shuffled_samples = copy.copy(sample_labels)
    cur_df = copy.copy(input_df).reset_index()

    for chrom in chroms:
        np.random.shuffle(shuffled_samples)
        renaming = {label: shuffled_label for label,
                    shuffled_label in zip(sample_labels, shuffled_samples)}
        renaming.update({normal_name: normal_name})
        cur_df.loc[cur_df['chrom'] == chrom, 'sample_id'] = cur_df.loc[cur_df['chrom']
                                                                       == chrom, 'sample_id'].map(lambda x: renaming[x])

    cur_df.sort_values(['sample_id', 'chrom', 'start', 'end'], inplace=True)
    cur_df.set_index(['sample_id', 'chrom', 'start', 'end'], inplace=True)
    cur_df.sort_index(inplace=True)

    return cur_df


def _bitstrs(tree):
    """Create a signature for the tree topology. Each internal node has a binary string with length
    equal to the number of leaf nodes with values 1 if the internal node is downstream of a particular
    leaf node and 0 otherwise.
    The final list of strings is a unique identifier of the tree topology.
    """
    # from https://biopython.org/wiki/Phylo
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString("".join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs


def compare_trees(tree1, tree2, fail_on_different_terminals=True):
    """Compare two tree topologies. The length of the branches are not regarded.

    If fail_on_different_terminals is set to True the function fails if the terminals are not the
    same. This prevents accidental comparison of the wrong trees.
    
    Modified from https://biopython.org/wiki/Phylo
    """
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        if fail_on_different_terminals:
            raise ValueError("Leaf names of the two trees are not the same:\ntree 1: {}\ntree 2: {}.".format(
                set(term_names1), set(term_names2)))
        else:
            return False
    # true if _BitStrings are the same
    return _bitstrs(tree1) == _bitstrs(tree2)


def _single_bootstrap_run(input_df, fst, bootstrap_method, i, N_bootstrap, normal_name='diploid'):
    cur_df = bootstrap_method(input_df)
    _, _, _, cur_final_tree, _, _ = main(
        input_df=cur_df,
        asymm_upper_fst=fst,
        asymm_lower_fst=fst,
        normal_name=normal_name,
        input_tree=None,
        ancestral_reconstruction=False,
        chr_separator='X')

    if (i+1) % (N_bootstrap//5) == 0:
        logger.info('{}/{} ({}%) bootstrap runs completed'.format(i + 1, N_bootstrap, 
                                                                  int((i+1)/N_bootstrap*100)))

    return cur_final_tree


def run_bootstrap(input_df,
                  original_tree=None,
                  N_bootstrap=50,
                  method='chr-wise',
                  wgd=True,
                  seed=None,
                  normal_name='diploid',
                  show_progress=True,
                  n_cores=None):
    """Run a given number of bootstrapping steps on the original data. 

    From the original data either a set of chromosome-wise bootstrap or segment-wise jackknife datasets
    are created. For every one of these datasets the optimal MEDICC tree is reconstructed and added
    to a tree dataframe. If the original tree is provided support values for the individual branches
    are calculated.
    """

    # tqdm can be used for progress bars
    try:
        from tqdm.auto import tqdm
    except ImportError:
        def tqdm(x, disable):
            return x

    if not show_progress:
        def tqdm(x, disable):
            return x

    if method == 'chr-wise':
        bootstrap_method = chr_wise_bootstrap_df
        logger.info('Starting {} chr-wise bootstrap runs'.format(N_bootstrap))
    elif method == 'segment-wise':
        logger.info('Starting {} segment-wise jackknife runs'.format(N_bootstrap))
        bootstrap_method = segment_wise_jacknife_df
    else:
        raise ValueError('method has to be either chr-wise or segment-wise.')

    if seed is not None:
        np.random.seed(seed)

    if wgd:
        fst = io.read_fst()
    else:
        fst = io.read_fst(no_wgd=True)

    if original_tree is not None:
        trees = {original_tree: [0, 'original']}
    else:
        trees = dict()

    # Run the actual bootstrapping steps
    if n_cores is not None and n_cores > 1:
        try:
            from joblib import Parallel, delayed
        except ImportError:
            raise ImportError("joblib must be installed for parallelization")
        initial_trees = Parallel(n_jobs=n_cores)(delayed(_single_bootstrap_run)(
            input_df, fst, bootstrap_method, i, N_bootstrap, normal_name)
            for i in range(N_bootstrap))
    else:
        initial_trees = []
        for i in tqdm(range(N_bootstrap), disable=not show_progress):
            cur_tree = _single_bootstrap_run(
                input_df, fst, bootstrap_method, i, N_bootstrap, normal_name)
            initial_trees.append(cur_tree)

    # delete duplicate trees
    for new_tree in initial_trees:
        for tree in trees.keys():
            if compare_trees(new_tree, tree):
                trees[tree][0] += 1
                break
        else:
            trees[new_tree] = [1, '']

    trees_df = pd.DataFrame(trees).T.reset_index().sort_values(
        0, ascending=False).reset_index(drop=True)
    trees_df.rename({0: 'count', 1: 'note', 'index': 'tree'}, axis=1, inplace=True)
    trees_df['freq'] = trees_df['count'] / trees_df['count'].sum()
    logger.info('A total of {} bootstrap trees were created'.format(len(trees_df)))

    # Calculate the support tree from the bootstrap trees
    if original_tree is not None:
        logger.info('Input tree provided. Calculating support tree')
        trees_list = []
        for _, row in trees_df.iterrows():
            trees_list += row['count'] * [row['tree']]
        support_tree = get_support(original_tree, trees_list)

        # Remove branch support for first branch because it is meaningless for rooted trees
        for clade in support_tree.root.clades:
            clade.confidence = None

    else:
        logger.debug('No input tree provided. Support tree is not calculated')
        support_tree = None

    return trees_df, support_tree


def bootstrap_hypothesis_test(input_df, hypothesis, trees_df=None, **kwargs):
    """Test a hypothesis on a set of bootstrap trees. `hypothesis` has to be a function.
    """
    if trees_df is None:
        trees_df, _ = run_bootstrap(input_df, **kwargs)

    trees_df['hypothesis'] = trees_df['tree'].map(hypothesis)

    return trees_df, trees_df.groupby('hypothesis')['freq'].sum()
