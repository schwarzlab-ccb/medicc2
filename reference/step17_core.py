# REFERENCE for migration Step 17 (not part of REF; untracked). Goes to medicc/core.py in NEW.
# It is NEW's committed file (master) with only this step's changes; see the plan for the step.
# main(): the fstengine path (strings, fstengine.pairwise_distances, reconstruct_ancestors(n_cores),
# create_df_from_cn_str_dict); --events raises MEDICCError until Phase 6 (Step 34).
# Helpers from REF (native -> fstengine): _check_fstengine_options, create_cn_str_dict,
# create_standard_fsa_dict_from_data (now uses create_cn_str_dict), create_df_from_fsa +
# create_df_from_cn_str_dict, update_branch_lengths (fstengine directed scores).
# detect_wgd() is deleted (as in REF): it used the removed read_fst(n_wgd=...) options.
# Not yet (later phases): FastME/tree_method (Steps 38-40), NNI + reconstructed_tree_length
# (Steps 36-37), engine event histories (Step 34). main() keeps chr_separator and wgd_x2 so
# bootstrap.py and medicc2 keep working until those steps remove them.

import copy
import logging
import os
from functools import lru_cache
from itertools import combinations

import Bio
import fstlib
import numpy as np
import pandas as pd

import medicc
from medicc import fstengine, io, nj, tools, event_reconstruction


# prepare logger 
logger = logging.getLogger(__name__)


def main(input_df,
         asymm_fst,
         normal_name='diploid',
         input_tree=None,
         ancestral_reconstruction=True,
         chr_separator='X',
         prune_weight=0,
         allele_columns=['cn_a', 'cn_b'],
         wgd_x2=False,
         no_wgd=False,
         total_cn=False,
         n_cores=None,
         reconstruct_events=False):
    """ MEDICC Main Method

    asymm_fst is the model from medicc.io.read_fst: the fstengine (medicc.fstengine.ModelConfig)
    or a custom fstlib FST. The fstengine requires prune_weight 0, and no_wgd=True restricts it
    to gains and losses. Event reconstruction (reconstruct_events=True) is not available yet:
    for the fstengine it comes with its event histories, a custom FST has no event labels.
    Raises MEDICCError for unsupported options.
    """
    use_fstengine = isinstance(asymm_fst, fstengine.ModelConfig)
    if reconstruct_events:
        if use_fstengine:
            raise MEDICCError("Event reconstruction is not available for the built-in model yet.")
        raise MEDICCError("Event reconstruction is not available for a custom FST "
                          "(an arbitrary FST has no event labels).")
    sample_labels = input_df.index.get_level_values('sample_id').unique()

    ## Validate input
    logger.info("Validating input.")
    if use_fstengine:
        asymm_fst = _check_fstengine_options(input_df, asymm_fst, prune_weight, no_wgd, total_cn)
        io.validate_input(input_df, normal_name=normal_name, min_samples=2)
        CN_str_dict = create_cn_str_dict(input_df)
        FSA_dict = CN_str_dict  # the fstengine works on copy-number strings, not FSAs
    else:
        symbol_table = asymm_fst.input_symbols()
        io.validate_input(input_df, symbol_table, normal_name=normal_name)

        ## Compile input data into FSAs stored in dictionaries
        logger.info("Compiling input sequences into FSAs.")
        FSA_dict, CN_str_dict = create_standard_fsa_dict_from_data(input_df, symbol_table, chr_separator)

    ## Reconstruct a tree
    if input_tree is None:
        ## Calculate pairwise distances
        logger.info("Calculating pairwise distance matrices")
        if use_fstengine:
            pairwise_distances = fstengine.pairwise_distances(CN_str_dict, asymm_fst, n_cores)
        elif n_cores is not None and n_cores > 1:
            pairwise_distances = parallelization_calc_pairwise_distance(sample_labels, asymm_fst, CN_str_dict,
                                                                                    n_cores)
        else:
            pairwise_distances = calc_pairwise_distance_matrix(asymm_fst, CN_str_dict)

        if (pairwise_distances == np.inf).any().any():
            affected_pairs = [(pairwise_distances.index[s1], pairwise_distances.index[s2])
                              for s1, s2 in zip(*np.where((pairwise_distances == np.inf)))]
            raise MEDICCError("Evolutionary distances could not be calculated for some sample "
                              "pairings. Please check the input data.\n\nThe affected pairs are: "
                              f"{affected_pairs}")

        logger.info("Inferring tree topology.")
        nj_tree = infer_tree_topology(
            pairwise_distances.values, pairwise_distances.index, normal_name=normal_name)
    else:
        logger.info("Tree provided, using it. No pairwise distance matrix is calculated!")

        pairwise_distances = pd.DataFrame(0, columns=FSA_dict.keys(), index=FSA_dict.keys())

        assert len([x for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name]) == \
            len(np.unique(input_df.index.get_level_values('sample_id'))), \
                "Number of samples differs in input tree and input dataframe"
        assert np.all(
            np.sort([x.name for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name]) ==
            np.sort(np.unique(input_df.index.get_level_values('sample_id')))), (
                "Input tree does not match input dataframe: "
                f"{np.sort([x.name for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name])}\n"
                f"{np.sort(np.unique(input_df.index.get_level_values('sample_id')))}")
        
        # necessary for the way that reconstruct_ancestors is performed
        if ancestral_reconstruction:
            input_tree.root_with_outgroup([x for x in input_tree.root.clades if x.name != normal_name][0].name)

        nj_tree = input_tree

    final_tree = copy.deepcopy(nj_tree)

    if ancestral_reconstruction:
        logger.info("Reconstructing ancestors.")
        ancestors = medicc.reconstruct_ancestors(tree=final_tree,
                                                 samples_dict=FSA_dict,
                                                 fst=asymm_fst,
                                                 normal_name=normal_name,
                                                 prune_weight=prune_weight,
                                                 n_cores=n_cores)

        ## Create and write output data frame with ancestors
        logger.info("Creating output copynumbers.")
        if use_fstengine:
            output_df = create_df_from_cn_str_dict(input_df, ancestors)
        else:
            output_df = create_df_from_fsa(input_df, ancestors)

        ## Update branch lengths with ancestors
        logger.info("Updating branch lengths of final tree using ancestors.")
        update_branch_lengths(final_tree, asymm_fst, ancestors, normal_name)
    else:
        output_df = input_df.copy()

    nj_tree.root_with_outgroup(normal_name)
    final_tree.root_with_outgroup(normal_name)

    if ancestral_reconstruction and reconstruct_events:
        logger.info("Reconstructing events.")
        output_df, events_df = event_reconstruction.calculate_all_cn_events(
            final_tree, output_df, allele_columns, normal_name,
            wgd_x2=wgd_x2, no_wgd=no_wgd, total_cn=total_cn)
        if len(events_df) != final_tree.total_branch_length():
            faulty_nodes = []
            for node in final_tree.find_clades():
                if node.name is not None and node.name != normal_name and node.branch_length != 0 and node.branch_length != len(events_df.loc[node.name]):
                    faulty_nodes.append(node.name)
            logger.warning("Event recreation was faulty. Events in '_cn_events_df.tsv' will be "
                        f"incorrect for the following nodes: {faulty_nodes}. "
                        f"total_branch_length: {final_tree.total_branch_length()}, "
                        f"nr of inferred events: {len(events_df)}")

    else:
        events_df = None

    return sample_labels, pairwise_distances, nj_tree, final_tree, output_df, events_df


def _check_fstengine_options(input_df, model_config, prune_weight, no_wgd, total_cn):
    """Check the main() options for the fstengine; returns the model config to use."""
    if prune_weight < 0:
        raise MEDICCError("prune_weight must not be negative.")
    if prune_weight > 0:
        raise MEDICCError("prune_weight > 0 is not supported by the built-in model yet (work in "
                          "progress). Use prune_weight=0.")
    if (total_cn or model_config.total_copy_numbers) and len(input_df.columns) != 1:
        raise MEDICCError("Total copy numbers require exactly one copy-number column.")
    if not set(np.unique(input_df.values)) <= set('012345678'):
        raise MEDICCError("The built-in model requires integer copy numbers from 0 to 8.")
    if no_wgd and not model_config.no_wgd:
        model_config = fstengine.ModelConfig(total_copy_numbers=model_config.total_copy_numbers,
                                             no_wgd=True)
    return model_config


def create_cn_str_dict(input_data, separator="X"):
    """ Creates a dictionary of copy-number strings from input DataFrame or Series.
    The keys of the dictionary are the sample/taxon names.
    If the input is a DataFrame, the string is the concatenated copy number profile of all allele
    columns, chromosomes separated by separator. """

    cn_str_dict = {}
    if isinstance(input_data, pd.DataFrame):
        logger.info('Creating copy-number strings for pd.DataFrame with the following data columns: {}'.format(
            input_data.columns.values))
        def aggregate_copy_number_profile(cnp):
            return separator.join([separator.join(["".join(x.astype('str'))
                                                   for _, x in cnp[allele].groupby('chrom', observed=False)]) for allele in cnp.columns])

    elif isinstance(input_data, pd.Series):
        logger.info('Creating copy-number strings for pd.Series with the name {}'.format(input_data.name))
        def aggregate_copy_number_profile(cnp):
            return separator.join(["".join(x.astype('str')) for _, x in cnp.groupby('chrom', observed=False)])

    else:
        raise MEDICCError("Input to function create_standard_fsa_dict_from_data has to be either"
                          "pd.DataFrame or pd.Series. \n input provided was {}".format(type(input_data)))
    
    for taxon, cnp in input_data.groupby('sample_id'):
        cn_str_dict[taxon] = aggregate_copy_number_profile(cnp)

    return cn_str_dict


def create_standard_fsa_dict_from_data(input_data,
                                       symbol_table: fstlib.SymbolTable,
                                       separator: str = "X") -> dict:
    """ Creates a dictionary of FSAs from input DataFrame or Series.
    The keys of the dictionary are the sample/taxon names. 
    If the input is a DataFrame, the FSA will be the concatenated copy number profiles of all allele columns"""

    cn_str_dict = create_cn_str_dict(input_data, separator)
    fsa_dict = {taxon: fstlib.factory.from_string(cn_str,
                                                  arc_type="standard",
                                                  isymbols=symbol_table,
                                                  osymbols=symbol_table)
                for taxon, cn_str in cn_str_dict.items()}

    return fsa_dict, cn_str_dict


def create_phasing_fsa_dict_from_df(input_df: pd.DataFrame, symbol_table: fstlib.SymbolTable, separator: str = "X") -> dict:
    """ Creates a dictionary of FSAs from two allele columns (Pandas DataFrame).
    The keys of the dictionary are the sample/taxon names. """
    allele_columns = input_df.columns
    if len(allele_columns) != 2:
        raise MEDICCError("Need exactly two alleles for phasing.")

    fsa_dict = {}
    for taxon, cnp in input_df.groupby('sample_id'):
        allele_a = cnp[allele_columns[0]]
        allele_b = cnp[allele_columns[1]]
        cn_str_a = separator.join(["".join(x) for _,x in allele_a.groupby(level='chrom', sort=False)])
        cn_str_b = separator.join(["".join(x) for _,x in allele_b.groupby(level='chrom', sort=False)])
        encoded = np.array([list(zip(cn_str_a, cn_str_b)), list(zip(cn_str_b, cn_str_a))])
        fsa_dict[taxon] = fstlib.factory.from_array(encoded, symbols=symbol_table, arc_type='standard')
        fsa_dict[taxon] = fstlib.determinize(fsa_dict[taxon]).minimize()

    return fsa_dict

def phase(input_df: pd.DataFrame, model_fst: fstlib.Fst, reference_sample='diploid', separator: str = 'X') -> pd.DataFrame:
    """ Phases every FST against the reference sample. 
    Returns two standard FSA dicts, one for each allele. """
    
    diploid_fsa = medicc.tools.create_diploid_fsa(model_fst)
    phasing_dict = medicc.create_phasing_fsa_dict_from_df(input_df, model_fst.input_symbols(), separator)
    fsa_dict_a, fsa_dict_b, _ = phase_dict(phasing_dict, model_fst, diploid_fsa)
    output_df = medicc.create_df_from_phasing_fsa(input_df, [fsa_dict_a, fsa_dict_b], separator)

    # Phasing across chromosomes is random, so we need to swap haplotype assignment per chromosome
    # so that the higher ploidy haplotype is always cn_a
    output_df['width'] = output_df.eval('end+1-start')
    output_df['cn_a_width'] = output_df['cn_a'].astype(float) * output_df['width']
    output_df['cn_b_width'] = output_df['cn_b'].astype(float) * output_df['width']

    swap_haplotypes_ind = output_df.groupby(['sample_id', 'chrom'])[
        ['cn_a_width', 'cn_b_width']].mean().diff(axis=1).iloc[:, 1] > 0

    output_df = output_df.join(swap_haplotypes_ind.rename('swap_haplotypes_ind'), on=['sample_id', 'chrom'])
    output_df.loc[output_df['swap_haplotypes_ind'], ['cn_a', 'cn_b']] = output_df.loc[output_df['swap_haplotypes_ind'], ['cn_b', 'cn_a']].values
    output_df = output_df.drop(['width', 'cn_a_width', 'cn_b_width', 'swap_haplotypes_ind'], axis=1)

    return output_df

def phase_dict(phasing_dict, model_fst, reference_fst):
    """ Phases every FST against the reference sample. 
    Returns two standard FSA dicts, one for each allele. """
    fsa_dict_a = {}    
    fsa_dict_b = {}
    scores = {}
    left = (reference_fst * model_fst).project('output')
    right = (~model_fst * reference_fst).project('input')
    for sample_id, sample_fst in phasing_dict.items():
        phased_fst = fstlib.align(sample_fst, left, right).topsort()
        score = fstlib.shortestdistance(phased_fst, reverse=True)[phased_fst.start()]
        scores[sample_id] = float(score)
        fsa_dict_a[sample_id] = fstlib.arcmap(phased_fst.copy().project('input'), map_type='rmweight')
        fsa_dict_b[sample_id] = fstlib.arcmap(phased_fst.project('output'), map_type='rmweight')
    
    return fsa_dict_a, fsa_dict_b, scores


def create_df_from_fsa(input_df: pd.DataFrame, fsa, separator: str = 'X'):
    """ 
    Takes a single FSA dict or a list of FSA dicts and extracts the copy number profiles.
    The allele names are taken from the input_df columns and the returned data frame has the same 
    number of rows and row index as the input_df. """

    if not isinstance(fsa, dict):
        raise MEDICCError("fsa input to create_df_from_fsa has to be a dict"
                          "Input type is {}".format(type(fsa)))

    samples = input_df.index.get_level_values('sample_id').unique()
    cn_str_dict = {node: tools.fsa_to_string(fsa[node]) for node in fsa if node not in samples}

    return create_df_from_cn_str_dict(input_df, cn_str_dict, separator)


def create_df_from_cn_str_dict(input_df: pd.DataFrame, cn_str_dict, separator: str = 'X'):
    """
    Takes a dict of copy-number strings (see create_cn_str_dict) and extracts the copy number
    profiles of the nodes that are not samples in input_df.
    The allele names are taken from the input_df columns and the returned data frame has the same
    number of rows and row index as the input_df. """

    alleles = input_df.columns
    nr_alleles = len(alleles)
    samples = input_df.index.get_level_values('sample_id').unique()
    output_df = input_df.unstack('sample_id')

    # Create dict and concat later to prevent pandas PerformanceWarning
    internal_cns = dict()
    for node in cn_str_dict:
        if node in samples:
            continue
        cns = cn_str_dict[node].split(separator)
        if len(cns) % nr_alleles != 0:
            raise MEDICCError('For sample {} we have {} haplotype-specific chromosomes for {} alleles'
                              '\nnumber of chromosomes has to be divisible by nr of alleles'.format(node,
                                                                                                    len(cns),
                                                                                                    nr_alleles))
        nr_chroms = int(len(cns) // nr_alleles)
        for i, allele in enumerate(alleles):
            cn = list(''.join(cns[(i*nr_chroms):((i+1)*nr_chroms)]))
            internal_cns[(allele, node)] = cn

    if len(internal_cns) == 0:  # e.g. a single tumour sample: no internal nodes
        return input_df.copy()

    internal_cns_df = pd.DataFrame(internal_cns, index=output_df.index)
    internal_cns_df.columns.names = ['allele', 'sample_id']
    output_df = (pd.concat([output_df, internal_cns_df], axis=1)
                 .stack('sample_id')
                 .reorder_levels(['sample_id', 'chrom', 'start', 'end'])
                 .sort_index())

    return output_df


def create_df_from_phasing_fsa(input_df: pd.DataFrame, fsas, separator: str = 'X'):
    """ 
    Takes a two FSAs dicts from phasing and extracts the copy number profiles.
    The allele names are taken from the input_df columns and the returned data frame has the same 
    number of rows and row index as the input_df. """

    alleles = input_df.columns
    if len(fsas) != 2:
        raise MEDICCError("fsas has to be of length 2")
    if not all([isinstance(fsa, dict) for fsa in fsas]):
        raise MEDICCError("all fsas entries have to be dicts")
    if fsas[0].keys() != fsas[1].keys():
        raise MEDICCError("fsas keys have to be the same")


    output_df = input_df.copy()[[]]
    output_df[alleles] = ''

    for sample in fsas[0].keys():
        cns_a = tools.fsa_to_string(fsas[0][sample]).split(separator)
        cns_b = tools.fsa_to_string(fsas[1][sample]).split(separator)
        if len(cns_a) != len(cns_b):
            raise MEDICCError(f"length of alleles is not the same for sample {sample}")

        output_df.loc[sample, alleles[0]] = list(''.join(cns_a))
        output_df.loc[sample, alleles[1]] = list(''.join(cns_b))

    # output_df = output_df.stack('sample_id')
    # output_df = output_df.reorder_levels(['sample_id', 'chrom', 'start', 'end']).sort_index()
    
    return output_df


def shorten_cn_strings(string_1, string_2):
    '''
    Takes two strings string_1 and string_2 and removes entires that are consecutive duplicates in both strings.

    Example:
        Input:
            string_1 = "abccd"
            string_2 = "1233d"
        Output:
            string_1_short = "abcd"
            string_2_short = "123d"
    '''
    assert len(string_1) == len(string_2)
    keep_indices = [i for i in range(len(string_1)) if
                    i == 0 or (string_1[i] != string_1[i - 1]) or (string_2[i] != string_2[i - 1])]
    string_1_short = ''.join(string_1[i] for i in keep_indices)
    string_2_short = ''.join(string_2[i] for i in keep_indices)
    return string_1_short, string_2_short


def parallelization_calc_pairwise_distance(sample_labels, asymm_fst, CN_str_dict, n_cores):
    try:
        from joblib import Parallel, delayed
    except ImportError:
        raise ImportError("joblib must be installed for parallelization")

    parallelization_groups = medicc.tools.create_parallelization_groups(len(sample_labels))
    parallelization_groups = [sample_labels[group] for group in parallelization_groups]
    logger.info("Running {} parallel runs on {} cores".format(len(parallelization_groups), n_cores))

    parallel_pairwise_distances = Parallel(n_jobs=n_cores)(delayed(calc_pairwise_distance_matrix)(
        asymm_fst, {key: val for key, val in CN_str_dict.items() if key in cur_group}, True)
            for cur_group in parallelization_groups)

    pdm = medicc.tools.total_pdm_from_parallel_pdms(sample_labels, parallel_pairwise_distances)

    return pdm


@lru_cache(maxsize=None)
def calc_MED_distance(model_fst, profile_1, profile_2):
    '''
    Calculate the MED distance between two profiles represented as strings.
    '''

    profile_1_short, profile_2_short = shorten_cn_strings(profile_1, profile_2)

    # Convert shrunken string to fsa
    symbol_table = model_fst.input_symbols()
    profile_1_short_fsa = fstlib.factory.from_string(profile_1_short, isymbols=symbol_table, osymbols=symbol_table)
    profile_2_short_fsa = fstlib.factory.from_string(profile_2_short, isymbols=symbol_table, osymbols=symbol_table)

    # Calculate the MED distance
    distance = float(fstlib.kernel_score(model_fst, profile_1_short_fsa, profile_2_short_fsa))

    return distance


def calc_pairwise_distance_matrix(model_fst, cn_str_dict, parallel_run=True):
    samples = list(cn_str_dict.keys())
    pdm = pd.DataFrame(0, index=samples, columns=samples, dtype=float)
    combs = list(combinations(samples, 2))
    ncombs = len(combs)

    for i, (sample_a, sample_b) in enumerate(combs):
        cur_dist = calc_MED_distance(model_fst, cn_str_dict[sample_a], cn_str_dict[sample_b])
        pdm.loc[sample_a, sample_b] = cur_dist
        pdm.loc[sample_b, sample_a] = cur_dist

        if not parallel_run and (100 * (i + 1) / ncombs) % 10 == 0:  # log every 10%
            logger.info(f'{(i + 1) / ncombs * 100:.2f}')

    return pdm


def infer_tree_topology(pairwise_distances, labels, normal_name):
    if len(labels) > 2:
        tree = nj.NeighbourJoining(pairwise_distances, labels).tree

        tmpsearch = [c for c in tree.find_clades(name = normal_name)]
        normal_node = tmpsearch[0]
        root_path = tree.get_path(normal_node)[::-1]

        if len(root_path)>1:
            new_root = root_path[1]
            tree.root_with_outgroup(new_root)
    else:
        clade_ancestor = Bio.Phylo.PhyloXML.Clade(branch_length=0, name='internal_1')
        clade_ancestor.clades = [Bio.Phylo.PhyloXML.Clade(
            name=label, branch_length=0 if label == normal_name else 1) for label in labels]

        tree = Bio.Phylo.PhyloXML.Phylogeny(root=clade_ancestor)
        tree.root_with_outgroup(normal_name)

    return tree


def update_branch_lengths(tree, fst, ancestor_fsa, normal_name='diploid'):
    """ Updates the branch lengths in the tree using the internal nodes supplied in the FSA dict.
    For the fstengine (fst is a medicc.fstengine.ModelConfig) the dict holds copy-number strings
    and the branch lengths are fstengine directed scores; an unreachable child raises MEDICCError.
    """
    if not isinstance(ancestor_fsa, dict):
        raise MEDICCError("input ancestor_fsa to function update_branch_lengths has to be either a dict"
                          "provided type is {}".format(type(ancestor_fsa)))

    if isinstance(fst, fstengine.ModelConfig):
        model = fstengine.Model(fst)

        def _distance_to_child(fst, ancestor_fsa, sample_1, sample_2):
            distance = model.score(ancestor_fsa[sample_1], ancestor_fsa[sample_2], directed=True)
            if not np.isfinite(distance):
                raise MEDICCError(f"{sample_2} cannot be reached from {sample_1}.")
            return distance
    else:
        model = None

        def _distance_to_child(fst, ancestor_fsa, sample_1, sample_2):
            return float(fstlib.score(fst, ancestor_fsa[sample_1], ancestor_fsa[sample_2]))

    if len(ancestor_fsa) == 2:
        child_clade = [x for x in tree.find_clades() if x.name is not None and x.name != normal_name][0]
        child_clade.branch_length = _distance_to_child(fst, ancestor_fsa, normal_name, child_clade.name)

    for clade in tree.find_clades():
        if clade.name is None:
            continue
        children = clade.clades
        if len(children) != 0:
            for child in children:
                if child.name == normal_name:  # exception: evolution goes from diploid to internal node
                    logger.debug(f'Updating MRCA branch length from {child.name} to {clade.name}')
                    brs = _distance_to_child(fst, ancestor_fsa, child.name, clade.name)
                else:
                    logger.debug(f'Updating branch length from {clade.name} to {child.name}')
                    brs = _distance_to_child(fst, ancestor_fsa, clade.name, child.name)
                logger.debug(f'branch length: {brs}')
                child.branch_length = brs

    if model is not None:
        model.close()


def summarize_patient(tree, pdm, sample_labels, normal_name='diploid', events_df=None):
    """Calculate several summary values for the provided samples

    Args:
        tree (Bio.Phylo.Tree): Phylogenetic tree
        pdm (pandas.DataFrame): Pairwise distance matrix between the samples
        sample_labels (list): List of all samples
        normal_name (str, optional): Name of normal sample. Defaults to 'diploid'.
        events_df (pandas.DataFrame, optional): DataFrame containg all copy-number events. Defaults to None.

    Returns:
        pandas.DataFrame: Summary DataFrame
    """    
    branch_lengths = []
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                branch_lengths.append(child.branch_length)

    nsamples = len(sample_labels)
    tree_length = np.sum(branch_lengths) if len(branch_lengths) > 0 else None
    avg_branch_length = np.mean(branch_lengths) if len(branch_lengths) > 0 else None
    min_branch_length = np.min(branch_lengths) if len(branch_lengths) > 0 else None
    max_branch_length = np.max(branch_lengths) if len(branch_lengths) > 0 else None
    median_branch_length = np.median(branch_lengths) if len(branch_lengths) > 0 else None
    # p_star = stats.star_topology_test(pdm)
    # p_clock = stats.molecular_clock_test(pdm,
    #                                      np.flatnonzero(np.array(sample_labels) == normal_name)[0])
    if events_df is None:
        wgd_status = "unknown (run with --events flag to detect WGDs)"
    else:
        if "wgd" in events_df['type'].values:
            wgd_status = "WGD on branch " + \
                "and ".join(events_df.loc[events_df['type'] ==
                                          'wgd'].index.get_level_values('sample_id'))
        else:
            wgd_status = "no WGD"

    result = pd.Series({
        'nsamples': nsamples,
        'normal_name': normal_name,
        'tree_length': tree_length,
        'mean_branch_length': avg_branch_length,
        'median_branch_length': median_branch_length,
        'min_branch_length': min_branch_length,
        'max_branch_length': max_branch_length,
        # 'p_star': p_star,
        # 'p_clock': p_clock,
        'wgd_status': wgd_status,
    })
    
    return result


class MEDICCError(Exception):
    pass
