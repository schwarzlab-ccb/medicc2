import logging
import re

import Bio
import fstlib
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

def set_sequences_on_tree_from_df(tree: Bio.Phylo.BaseTree, df: pd.DataFrame, clear_before=True):
    """DEPRECATED!
    
    Set sequences on tree from dataframe

    Args:
        tree (Bio.Phylo.BaseTree): Tree to set sequences on
        df (pd.DataFrame): DataFrame with copy number information
        clear_before (bool, optional): Clear old sequences. Defaults to True.
    """
    raise DeprecationWarning("This function is deprecated and is not used anymore.")

    if not hasattr(tree.root, 'sequences'):
        tree = tree.as_phyloxml()
    for clade in tree.find_clades():
        if clear_before:
            clade.sequences.clear()
        for label, data in df.items():
            try:
                clade.sequences.append(
                    Bio.Phylo.PhyloXML.Sequence(
                        name='X'.join(data.loc[clade.name].groupby('chrom').apply(lambda x: ''.join(x))),
                        symbol=label.upper())) 
            except KeyError:
                pass


def fsa_to_string(fsa):
    fsa_string = fstlib.tools.strings(fsa)['input'].values[0]
    return fsa_string

def hex2int(x):
    return int(x, 16)

def int2hex(x):
    return hex(x)[2:]

def format_chromosomes(ds):
    """ Expects pandas Series with chromosome names. 
    The goal is to take recognisable chromosome names, i.e. chr4 or chrom3 and turn them into chr3 format.
    If the chromosomes names are not recognized, return them unchanged."""

    ds = ds.astype('str')
    pattern = re.compile(r"(chr|chrom)?(_)?(0)?((\d+)|X|Y)", flags=re.IGNORECASE)
    matches = ds.apply(pattern.match)
    matchable = ~matches.isnull().any()
    if matchable:
        newchr = matches.apply(lambda x:f"chr{x[4].upper():s}")
        numchr = matches.apply(lambda x:int(x[5]) if x[5] is not None else -1)
        chrlevels = np.sort(numchr.unique())
        chrlevels = np.setdiff1d(chrlevels, [-1])
        chrcats = [f"chr{i}" for i in chrlevels]
        if 'chrX' in list(newchr):
            chrcats += ['chrX',]
        if 'chrY' in list(newchr):
            chrcats += ['chrY',]
        newchr = pd.Categorical(newchr, categories=chrcats)
    else:
        logger.warning("Could not match the chromosome labels. Rename the chromosomes to follow the chr1, "
                    "chr2, ... convention to avoid potential errors. "
                    "Current format: {}.".format(ds.unique()))
        newchr = pd.Categorical(ds, categories=ds.unique())
    assert not newchr.isna().any(), "Could not reformat chromosome labels. Rename according to chr1, chr2, ..."
    return newchr


def next_prime(N):
    '''calculate a prime number p with p >= N'''
    def is_prime(x):
        return all(x % i for i in range(3, x))

    cur_p = int(np.ceil(np.sqrt(N)))
    if cur_p % 2 == 0:
        cur_p += 1
    while not is_prime(cur_p):
        cur_p += 2
    return cur_p


def create_parallelization_groups(number_samples):
    '''Create subgroups of the samples to perform parallel runs of the pairwise distances matrix 
    calculations. 

    Args:
        number_samples (int): Total number of samples

    Returns:
        list of lists: Lists, each one containing up to p indices


    Method implemented as proposed in:
    Emmanuel Sapin, Matthew C Keller, Novel approach for parallelizing pairwise comparison problems as applied to detecting segments identical by descent in whole-genome data,
    Bioinformatics, 2021;, btab084, https://doi.org/10.1093/bioinformatics/btab084
    '''
    p = next_prime(number_samples)
    p_matrix = np.arange(p**2).reshape((p, p)).T
    p_matrix[p_matrix >= number_samples] = -1

    p_matrix_permutations = [p_matrix]
    for _ in range(p-1):
        p_matrix_permutations.append(np.array([np.roll(p_matrix_permutations[-1][:, i], -i) for i in range(p)]).T)

    groups = [*p_matrix.T] + [x for cur_p in p_matrix_permutations for x in [*cur_p]]
    # remove -1 entries (numbers larger than number_samples)
    groups = [group[group >= 0] for group in groups]
    groups = [group for group in groups if len(group) > 1]

    return groups


def total_pdm_from_parallel_pdms(sample_labels, parallel_pdms):
    total_pdm = pd.DataFrame(index=sample_labels, columns=sample_labels, dtype=float)

    for cur_pdm in parallel_pdms:
        total_pdm.loc[cur_pdm.index, cur_pdm.index] = cur_pdm

    if total_pdm.isna().sum().sum() != 0:
        raise ValueError('Could not compute pairwise distances for these sample pairs indices:\n{}'.format(
            np.where(total_pdm.isna())))

    return total_pdm


def create_diploid_fsa(fst, total_copy_numbers=False):

    # Create diploid FSA
    diploid_fsa = fstlib.Fst()
    diploid_fsa.set_input_symbols(fst.input_symbols())
    diploid_fsa.set_output_symbols(fst.output_symbols())
    diploid_fsa.add_state()
    diploid_fsa.set_start(0)
    diploid_fsa.set_final(0, 0)
    if total_copy_numbers:
        diploid_fsa.add_arc(0, ('2', '2', 0, 0))
    else:
        diploid_fsa.add_arc(0, ('1', '1', 0, 0))
    diploid_fsa.add_arc(0, ('X', 'X', 0, 0))

    return diploid_fsa

def sum_of_branch_length(tree):
    '''Calculate the sum of all branch lengths in the tree'''
    return sum([clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None])