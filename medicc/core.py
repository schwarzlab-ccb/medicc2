import copy
import logging
import os
import random
from functools import lru_cache
from itertools import combinations

import Bio
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import fstlib
import numpy as np
import pandas as pd

import medicc
from medicc import io, nj, tools, event_reconstruction, tree_hash
from medicc.tools import int2hex

import matplotlib.pyplot as plt

# prepare logger
logger = logging.getLogger(__name__)


def main(input_df,
         asymm_fst,
         output_dir,
         event_counting_fst,
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
         reconstruct_events=False,
         euclidean=False,
         test_flag_for_length_encoding_only_on_the_lower_half_delete_when_release_to_public=False,
         constrained_sankoff=False):
    """ MEDICC Main Method """

    symbol_table = asymm_fst.input_symbols()

    ## Validate input
    logger.info("Validating input.")
    io.validate_input(input_df, symbol_table, normal_name=normal_name)

    ## Compile input data into FSAs stored in dictionaries
    logger.info("Compiling input sequences into FSAs.")
    FSA_dict, CN_str_dict = create_standard_fsa_dict_from_data(input_df, symbol_table, chr_separator)
    sample_labels = input_df.index.get_level_values('sample_id').unique()

    ## Reconstruct a tree
    if input_tree is None:
        if n_cores is not None and n_cores > 1:
            pairwise_distances = parallelization_calc_pairwise_distance(sample_labels, asymm_fst, CN_str_dict,
                                                                                    n_cores, euclidean = euclidean, chr_separator = chr_separator)
        else:
            pairwise_distances = calc_pairwise_distance_matrix(asymm_fst, CN_str_dict, euclidean=euclidean, chr_separator=chr_separator)

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

        tree_leaves = [x.name for x in list(input_tree.find_clades())
                       if x.name is not None and 'internal' not in x.name and x.name != normal_name]
        df_samples = np.unique(input_df.index.get_level_values('sample_id'))[1:]

        assert len(tree_leaves) == len(df_samples), \
            "Number of samples differs in input tree and input dataframe"

        assert np.all(np.sort(tree_leaves) == np.sort(df_samples)), (
            "Input tree does not match input dataframe: "
            f"{np.sort(tree_leaves)}\n"
            f"{np.sort(df_samples)}")
        
        # necessary for the way that reconstruct_ancestors is performed
        if ancestral_reconstruction:
            input_tree.root_with_outgroup([x for x in input_tree.root.clades if x.name != normal_name][0].name)

        nj_tree = input_tree

    final_tree = copy.deepcopy(nj_tree)

    if ancestral_reconstruction:
        logger.info("Reconstructing ancestors.")
        # default pass use the user/default fst for upper and lower pass
        upper_fst = asymm_fst
        lower_pass_fst = asymm_fst

        if test_flag_for_length_encoding_only_on_the_lower_half_delete_when_release_to_public:
            upper_fst = event_counting_fst
            lower_pass_fst = asymm_fst

        if constrained_sankoff:
            # event-minimal candidate sets (event FST) + length-minimising Sankoff
            # DP under the user-provided length-encoding FST (use a high open value,
            # e.g. open=5000, so the result stays event-minimal).
            upper_fst = event_counting_fst
            lower_pass_fst = asymm_fst
            ancestors, _uppass_cache = medicc.reconstruct_ancestors_constrained_sankoff(
                                                 tree=final_tree,
                                                 samples_dict=FSA_dict,
                                                 upper_pass_fst=upper_fst,
                                                 lower_pass_fst=lower_pass_fst,
                                                 normal_name=normal_name,
                                                 prune_weight=prune_weight,
                                                 spr_logger_disable=False,
                                                 n_cores=n_cores)
        else:
            ancestors, _uppass_cache = medicc.reconstruct_ancestors(tree=final_tree,
                                                 samples_dict=FSA_dict,
                                                 upper_pass_fst=upper_fst,
                                                 lower_pass_fst=lower_pass_fst,
                                                 normal_name=normal_name,
                                                 prune_weight=prune_weight,
                                                 spr_logger_disable=False,
                                                 n_cores=n_cores)

        ## Create and write output data frame with ancestors
        logger.info("Creating output copynumbers.")
        output_df = create_df_from_fsa(input_df, ancestors)

        ## Update branch lengths with ancestors
        logger.info("Updating branch lengths of final tree using ancestors.")
        update_branch_lengths(final_tree, event_counting_fst, ancestors, normal_name)
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


def _reshape_nj_for_search(nj_tree, normal_name):
    """Reshape an NJ tree for in-place tree search.

    After this transformation:
      - tree.root.name == normal_name
      - tree.root has exactly one structural child (the MRCA, e.g. "internal_0")
      - the diploid leaf is no longer a child of the root

    Mutates and returns the input tree.
    """
    nj_tree.root_with_outgroup(normal_name)
    nj_tree.root.clades = [clade for clade in nj_tree.root.clades if clade.name != normal_name]
    nj_tree.root.name = normal_name
    return nj_tree


def _wrap_tree_for_output(final_tree, fst, ancestors, normal_name):
    """Wrap a search-shape tree back into the standard PhyloXML output form.

    The search shape has tree.root named after the diploid with a single
    structural child (the MRCA). For plotting and output, MEDICC2 expects a new
    unnamed root with two clades: the diploid and the MRCA. This function does
    that wrap and updates branch lengths using the supplied FST and ancestors.

    Mutates and returns the input tree.
    """
    new_root_clade = Bio.Phylo.PhyloXML.Clade(branch_length=0)
    final_tree.root.branch_length = 0
    new_root_clade.clades.append(final_tree.root)
    new_root_clade.clades.append(final_tree.root.clades[0])
    final_tree.root.clades = []
    final_tree.root = new_root_clade

    logger.info("Updating branch lengths of final tree using ancestors.")
    update_branch_lengths(final_tree, fst, ancestors, normal_name)
    return final_tree


def main_spr(input_df,
             asymm_fst,
             output_dir,
             event_counting_fst,
             normal_name='diploid',
             input_tree=None,
             chr_separator='X',
             n_cores=None,
             prune_weight=0,
             spr_step=1000,
             spr_start="nj",
             sa_temp_start=10.0,
             sa_temp_end=0.01,
             sa_cooling="geometric"):
    """
    spr mode Main Method
    We start with a
        - if spr_start = "random" : A random tree
        - if spr_start = "neighbor-joining" : A neighbor-joining tree
        - if input_tree is not None : A user-provided tree, medicc3_MCMC_start is ignored in this case
    And uses SPR tree search moves to explore the tree space.
    """

    symbol_table = asymm_fst.input_symbols()

    ## Validate input
    logger.info("Validating input.")
    io.validate_input(input_df, symbol_table, normal_name=normal_name)

    ## Compile input data into FSAs stored in dictionaries
    logger.info("Compiling input sequences into FSAs.")
    FSA_dict, CN_str_dict = create_standard_fsa_dict_from_data(input_df, symbol_table, chr_separator)
    sample_labels = input_df.index.get_level_values('sample_id').unique()

    use_multichain = n_cores is not None and n_cores > 1
    n_chains = n_cores if use_multichain else 1
    chain_start_trees = []

    ## Tree Reconstruction

    ### Create start point tree
    if input_tree is None:
        if spr_start == "random":
            logger.info("SPR mode: Start with a random tree.")
            if use_multichain:
                logger.info("SPR mode: Multi-chain enabled with %s chains. Generating random start tree per chain.",
                            n_chains)
                chain_start_trees = _create_unique_random_start_trees(
                    sample_labels=sample_labels,
                    normal_name=normal_name,
                    n_chains=n_chains,
                )
                nj_tree = copy.deepcopy(chain_start_trees[0])
            else:
                nj_tree = create_random_tree_topology(list(sample_labels), normal_name=normal_name)
                chain_start_trees = [copy.deepcopy(nj_tree)]
        elif spr_start == "neighbor-joining" or spr_start == "nj":
            logger.info("SPR mode: Start with a Neighbor-joining tree.")
            logger.info("SPR mode: Calculating pairwise distance matrices.")
            if n_cores is not None and n_cores > 1:
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

            logger.info("SPR mode: Inferring tree topology using neighbor-joining.")
            nj_tree = infer_tree_topology(
                pairwise_distances.values, pairwise_distances.index, normal_name=normal_name)

            logger.debug(
                "SPR mode: Adjust the Neighbor-joining tree structure to be compatible with the MCMC process.")
            nj_tree = _reshape_nj_for_search(nj_tree, normal_name)

            if use_multichain:
                logger.info("SPR mode: Multi-chain enabled with %s chains. All chains start from the same neighbor-joining tree.",
                            n_chains)
                chain_start_trees = [copy.deepcopy(nj_tree) for _ in range(n_chains)]
            else:
                chain_start_trees = [copy.deepcopy(nj_tree)]
        else:
            raise NotImplementedError("SPR mode: Start with a random tree or a neighbor-joining tree.")
    else:
        logger.info("SPR mode: Tree provided, using it as the starting point.")

        assert len([x for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name]) == \
               len(np.unique(input_df.index.get_level_values('sample_id'))), \
            "Number of samples differs in input tree and input dataframe"
        assert np.all(
            np.sort(
                [x.name for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name]) ==
            np.sort(np.unique(input_df.index.get_level_values('sample_id')))), (
            "Input tree does not match input dataframe: "
            f"{np.sort([x.name for x in list(input_tree.find_clades()) if x.name is not None and 'internal' not in x.name])}\n"
            f"{np.sort(np.unique(input_df.index.get_level_values('sample_id')))}")

        input_tree.root_with_outgroup([x for x in input_tree.root.clades if x.name != normal_name][0].name)
        nj_tree = input_tree

        if use_multichain:
            logger.info("SPR mode: Multi-chain enabled with %s chains. All chains start from the provided input tree.",
                        n_chains)
            chain_start_trees = [copy.deepcopy(nj_tree) for _ in range(n_chains)]
        else:
            chain_start_trees = [copy.deepcopy(nj_tree)]

    logger.info("SPR mode: Running MCMC to infer phylogenetic tree")
    final_tree_l, ancestors_l, chain_traces = spr_mode(tree=copy.deepcopy(chain_start_trees[0]),
                                                       samples_dict=FSA_dict,
                                                       fst=asymm_fst,
                                                       normal_name=normal_name,
                                                       prune_weight=prune_weight,
                                                       MCMC_step=spr_step,
                                                       n_cores=n_cores,
                                                       sa_temp_start=sa_temp_start,
                                                       sa_temp_end=sa_temp_end,
                                                       sa_cooling=sa_cooling,
                                                       chain_start_trees=chain_start_trees)
    logger.info("SPR mode: Creating output copynumbers.")
    output_df_l = [create_df_from_fsa(input_df, ancestors) for ancestors in ancestors_l]

    # Adjust the tree format for MEDICC2's plotting function
    for i, final_tree in enumerate(final_tree_l):
        ancestors = ancestors_l[i]
        _wrap_tree_for_output(final_tree, event_counting_fst, ancestors, normal_name)




    # do the same to the nj_tree
    new_root_clade = Bio.Phylo.PhyloXML.Clade(branch_length=0)
    nj_tree.root.branch_length = 0
    new_root_clade.clades.append(nj_tree.root)
    new_root_clade.clades.append(nj_tree.root.clades[0])
    nj_tree.root.clades = []
    nj_tree.root = new_root_clade

    logger.info("SPR mode: Plotting SPR Trace")
    is_multichain_trace = len(chain_traces) > 0 and isinstance(chain_traces[0], list)
    chain_trace_l = chain_traces if is_multichain_trace else [chain_traces]

    # Save overall summary trace (single file).
    fig, ax = plt.subplots(figsize=(10, 6))
    if len(chain_trace_l) > 1:
        for chain_idx, chain_trace in enumerate(chain_trace_l, start=1):
            ax.plot(chain_trace, alpha=0.55, linewidth=1.2, label=f'chain_{chain_idx}')
        if len(chain_trace_l) <= 12:
            ax.legend(fontsize=8, ncol=2)
        ax.set_title('SPR Trace Summary (Multi-chain)')
    else:
        ax.plot(chain_trace_l[0], '-o', alpha=0.6)
        ax.set_title('SPR Trace Summary')
    ax.set_xlabel('SPR Step')
    ax.set_ylabel('Sum of Branch Lengths')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'SPR_trace.pdf'), bbox_inches='tight')
    plt.close()

    # Save one trajectory plot per chain.
    trace_output_dir = os.path.join(output_dir, 'spr_traces')
    os.makedirs(trace_output_dir, exist_ok=True)
    for chain_idx, chain_trace in enumerate(chain_trace_l, start=1):
        fig_chain, ax_chain = plt.subplots(figsize=(10, 6))
        ax_chain.plot(chain_trace, '-o', alpha=0.7, markersize=3)
        ax_chain.set_xlabel('SPR Step')
        ax_chain.set_ylabel('Sum of Branch Lengths')
        ax_chain.set_title(f'SPR Trace Chain {chain_idx}')
        ax_chain.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(os.path.join(trace_output_dir, f'SPR_trace_chain_{chain_idx}.pdf'),
                    bbox_inches='tight')
        plt.close()

    return sample_labels, nj_tree, final_tree_l, output_df_l


def main_nni(input_df,
             asymm_fst,
             output_dir,
             event_counting_fst,
             normal_name='diploid',
             input_tree=None,
             chr_separator='X',
             n_cores=None,
             prune_weight=0,
             nni_max_iter=100,
             nni_trace_dir=None,
             nni_start="nj"):
    """NNI hill-climbing main method.

    Starts from the MEDICC NJ tree (or a random tree if nni_start='random', or
    `input_tree` if provided) and runs deterministic steepest-ascent NNI
    hill-climbing until no neighbor strictly improves the score, or
    `nni_max_iter` sweeps are reached.

    nni_start: Starting tree for NNI search. One of 'nj' / 'neighbor-joining'
        (default) or 'random'. Ignored when input_tree is provided.
    nni_trace_dir: If not None, write one file per evaluated NNI neighbor into
        this directory, then aggregate into nni_trace.tsv.
    """
    symbol_table = asymm_fst.input_symbols()

    logger.info("Validating input.")
    io.validate_input(input_df, symbol_table, normal_name=normal_name)

    logger.info("Compiling input sequences into FSAs.")
    FSA_dict, CN_str_dict = create_standard_fsa_dict_from_data(input_df, symbol_table, chr_separator)
    sample_labels = input_df.index.get_level_values('sample_id').unique()

    if input_tree is None:
        if nni_start == "random":
            logger.info("NNI mode: Start with a random tree.")
            nj_tree = create_random_tree_topology(list(sample_labels), normal_name=normal_name)
        elif nni_start in ("nj", "neighbor-joining"):
            logger.info("NNI mode: Start with a Neighbor-joining tree.")
            logger.info("NNI mode: Calculating pairwise distance matrices.")
            if n_cores is not None and n_cores > 1:
                pairwise_distances = parallelization_calc_pairwise_distance(
                    sample_labels, asymm_fst, CN_str_dict, n_cores)
            else:
                pairwise_distances = calc_pairwise_distance_matrix(asymm_fst, CN_str_dict)

            if (pairwise_distances == np.inf).any().any():
                affected_pairs = [(pairwise_distances.index[s1], pairwise_distances.index[s2])
                                  for s1, s2 in zip(*np.where((pairwise_distances == np.inf)))]
                raise MEDICCError("Evolutionary distances could not be calculated for some sample "
                                  "pairings. Please check the input data.\n\nThe affected pairs are: "
                                  f"{affected_pairs}")

            logger.info("NNI mode: Inferring tree topology using neighbor-joining.")
            nj_tree = infer_tree_topology(
                pairwise_distances.values, pairwise_distances.index, normal_name=normal_name)
            logger.debug("NNI mode: Adjust the Neighbor-joining tree structure for tree search.")
            nj_tree = _reshape_nj_for_search(nj_tree, normal_name)
        else:
            raise NotImplementedError(f"NNI mode: Unknown start option '{nni_start}'. "
                                      "Use 'nj', 'neighbor-joining', or 'random'.")
    else:
        logger.info("NNI mode: Tree provided, using it as the starting point.")

        tree_leaves = [x.name for x in list(input_tree.find_clades())
                       if x.name is not None and 'internal' not in x.name and x.name != normal_name]
        df_samples = np.unique(input_df.index.get_level_values('sample_id'))[1:]

        assert len(tree_leaves) == len(df_samples), \
            "Number of samples differs in input tree and input dataframe"
        assert np.all(np.sort(tree_leaves) == np.sort(df_samples)), (
            "Input tree does not match input dataframe: "
            f"{np.sort(tree_leaves)}\n{np.sort(df_samples)}")

        input_tree.root_with_outgroup(
            [x for x in input_tree.root.clades if x.name != normal_name][0].name)
        nj_tree = input_tree

    logger.info("NNI mode: Running NNI hill-climbing to infer phylogenetic tree.")
    nni_result = nni_mode(
        tree=copy.deepcopy(nj_tree),
        samples_dict=FSA_dict,
        fst=asymm_fst,
        normal_name=normal_name,
        prune_weight=prune_weight,
        nni_max_iter=nni_max_iter,
        n_cores=n_cores,
        nni_trace_dir=nni_trace_dir,
    )
    final_tree_l = nni_result["best_trees"]
    ancestors_l = nni_result["best_ancestors"]
    trace = nni_result["trace"]
    step_records = nni_result["step_records"]

    logger.info("NNI mode: Creating output copynumbers.")
    output_df_l = [create_df_from_fsa(input_df, ancestors) for ancestors in ancestors_l]

    for i, final_tree in enumerate(final_tree_l):
        ancestors = ancestors_l[i]
        _wrap_tree_for_output(final_tree, event_counting_fst, ancestors, normal_name)

    # Wrap the nj_tree the same way (mirrors main_spr)
    new_root_clade = Bio.Phylo.PhyloXML.Clade(branch_length=0)
    nj_tree.root.branch_length = 0
    new_root_clade.clades.append(nj_tree.root)
    new_root_clade.clades.append(nj_tree.root.clades[0])
    nj_tree.root.clades = []
    nj_tree.root = new_root_clade

    logger.info("NNI mode: Plotting NNI trace.")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(trace, '-o', alpha=0.7, markersize=4)
    ax.set_xlabel('NNI Sweep')
    ax.set_ylabel('Sum of Branch Lengths')
    ax.set_title('NNI Trace')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'NNI_trace.pdf'), bbox_inches='tight')
    plt.close()

    if nni_trace_dir is not None and step_records:
        logger.info("NNI mode: Writing per-step trace files to %s", nni_trace_dir)
        os.makedirs(nni_trace_dir, exist_ok=True)
        # Per-step .txt files are the primary artifacts; nni_trace.tsv is aggregated from them for convenience.
        for step, newick_str, score in step_records:
            fname = os.path.join(nni_trace_dir, f"step_{step:08d}.txt")
            with open(fname, "w") as f:
                f.write(newick_str + "\n")
                f.write(str(float(score)) + "\n")

        logger.info("NNI mode: Aggregating trace files into nni_trace.tsv")
        step_files = sorted(
            f for f in os.listdir(nni_trace_dir) if f.startswith("step_") and f.endswith(".txt")
        )
        tsv_path = os.path.join(nni_trace_dir, "nni_trace.tsv")
        with open(tsv_path, "w") as tsv:
            tsv.write("step\tnewick\tsum_of_branch_length\n")
            for fname in step_files:
                step_num = int(fname[len("step_"):-len(".txt")])
                with open(os.path.join(nni_trace_dir, fname)) as f:
                    lines = f.read().strip().splitlines()
                newick_str = lines[0]
                score_str = lines[1]
                tsv.write(f"{step_num}\t{newick_str}\t{score_str}\n")

    return sample_labels, nj_tree, final_tree_l, output_df_l


def create_standard_fsa_dict_from_data(input_data,
                                       symbol_table: fstlib.SymbolTable,
                                       separator: str = "X") -> dict:
    """ Creates a dictionary of FSAs from input DataFrame or Series.
    The keys of the dictionary are the sample/taxon names. 
    If the input is a DataFrame, the FSA will be the concatenated copy number profiles of all allele columns"""

    fsa_dict = {}
    cn_str_dict = {}
    if isinstance(input_data, pd.DataFrame):
        logger.info('Creating FSA for pd.DataFrame with the following data columns: {}'.format(
            input_data.columns.values))
        def aggregate_copy_number_profile(cnp):
            return separator.join([separator.join(["".join(x.astype('str'))
                                                   for _, x in cnp[allele].groupby('chrom', observed=False)]) for allele in cnp.columns])

    elif isinstance(input_data, pd.Series):
        logger.info('Creating FSA for pd.Series with the name {}'.format(input_data.name))
        def aggregate_copy_number_profile(cnp):
            return separator.join(["".join(x.astype('str')) for _, x in cnp.groupby('chrom', observed=False)])

    else:
        raise MEDICCError("Input to function create_standard_fsa_dict_from_data has to be either"
                          "pd.DataFrame or pd.Series. \n input provided was {}".format(type(input_data)))
    
    for taxon, cnp in input_data.groupby('sample_id'):
        cn_str = aggregate_copy_number_profile(cnp)
        fsa_dict[taxon] = fstlib.factory.from_string(cn_str,
                                                     arc_type="standard",
                                                     isymbols=symbol_table,
                                                     osymbols=symbol_table)
        cn_str_dict[taxon] = cn_str

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

    alleles = input_df.columns
    if not isinstance(fsa, dict):
        raise MEDICCError("fsa input to create_df_from_fsa has to be a dict"
                          "Input type is {}".format(type(fsa)))

    nr_alleles = len(alleles)
    samples = input_df.index.get_level_values('sample_id').unique()
    output_df = input_df.unstack('sample_id')

    # Create dict and concat later to prevent pandas PerformanceWarning
    internal_cns = dict()
    for node in fsa:
        if node in samples:
            continue
        cns = tools.fsa_to_string(fsa[node]).split(separator)
        if len(cns) % nr_alleles != 0:
            raise MEDICCError('For sample {} we have {} haplotype-specific chromosomes for {} alleles'
                              '\nnumber of chromosomes has to be divisible by nr of alleles'.format(node,
                                                                                                    len(cns),
                                                                                                    nr_alleles))
        nr_chroms = int(len(cns) // nr_alleles)
        for i, allele in enumerate(alleles):
            cn = list(''.join(cns[(i*nr_chroms):((i+1)*nr_chroms)]))
            internal_cns[(allele, node)] = cn

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


def parallelization_calc_pairwise_distance(sample_labels, asymm_fst, CN_str_dict, n_cores, euclidean=False, chr_separator="X"):
    try:
        from joblib import Parallel, delayed
    except ImportError:
        raise ImportError("joblib must be installed for parallelization")

    parallelization_groups = medicc.tools.create_parallelization_groups(len(sample_labels))
    parallelization_groups = [sample_labels[group] for group in parallelization_groups]
    logger.info("Running {} parallel runs on {} cores".format(len(parallelization_groups), n_cores))

    parallel_pairwise_distances = Parallel(n_jobs=n_cores)(
        delayed(calc_pairwise_distance_matrix)(
            asymm_fst,
            {key: val for key, val in CN_str_dict.items() if key in cur_group},
            True,
            euclidean,
            chr_separator
        ) for cur_group in parallelization_groups)

    pdm = medicc.tools.total_pdm_from_parallel_pdms(sample_labels, parallel_pairwise_distances)

    return pdm


@lru_cache(maxsize=None)
def calc_MED_distance(model_fst, profile_1, profile_2, chr_separator="X", euclidean=False, disable_shortening=True):
    '''
    Calculate the MED distance between two profiles represented as strings.
    '''

    if not euclidean:
        profile_1_list = profile_1.split(chr_separator)
        profile_2_list = profile_2.split(chr_separator)
        profile_1_list_wrap_each_chr_with_sep = [chr_separator + chrom + chr_separator for chrom in profile_1_list]
        profile_2_list_wrap_each_chr_with_sep = [chr_separator + chrom + chr_separator for chrom in profile_2_list]
        profile_1_wrapped = ''.join(profile_1_list_wrap_each_chr_with_sep)
        profile_2_wrapped = ''.join(profile_2_list_wrap_each_chr_with_sep)
        profile_1_short = profile_1_wrapped
        profile_2_short = profile_2_wrapped
        if not disable_shortening:
            profile_1_short, profile_2_short = shorten_cn_strings(profile_1, profile_2)
        # Convert shrunken string to fsa
        symbol_table = model_fst.input_symbols()
        profile_1_short_fsa = fstlib.factory.from_string(profile_1_short, isymbols=symbol_table, osymbols=symbol_table)
        profile_2_short_fsa = fstlib.factory.from_string(profile_2_short, isymbols=symbol_table, osymbols=symbol_table)
        # Calculate the MED distance
        distance = float(fstlib.kernel_score(model_fst, profile_1_short_fsa, profile_2_short_fsa))
    else:
        profile_1_no_sep = profile_1.replace(chr_separator, '')
        profile_2_no_sep = profile_2.replace(chr_separator, '')

        profile_1_np_array = np.array([int(x) for x in profile_1_no_sep])
        profile_2_np_array = np.array([int(x) for x in profile_2_no_sep])

        distance = np.linalg.norm(profile_1_np_array - profile_2_np_array)

    return distance


def calc_pairwise_distance_matrix(model_fst, cn_str_dict, parallel_run=True, euclidean=False, chr_separator="X"):
    samples = list(cn_str_dict.keys())
    pdm = pd.DataFrame(0, index=samples, columns=samples, dtype=float)
    combs = list(combinations(samples, 2))
    ncombs = len(combs)

    for i, (sample_a, sample_b) in enumerate(combs):
        cur_dist = calc_MED_distance(model_fst, cn_str_dict[sample_a], cn_str_dict[sample_b], euclidean = euclidean, chr_separator = chr_separator)
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

def create_random_tree_topology(sample_labels, normal_name="diploid"):
    '''Create a random tree topology with the given sample labels and a normal sample name.'''

    # remove normal_name from sample_labels
    if normal_name in sample_labels:
        sample_labels.remove(normal_name)

    # Create terminal clades for each taxon
    clades = [Clade(name=label) for label in sample_labels]

    internal_counter = 1
    while len(clades) > 1:
        # Pick two random clades
        i1 = np.random.choice(range(len(clades)))
        clade1 = clades.pop(i1)
        i2 = np.random.choice(range(len(clades)))
        clade2 = clades.pop(i2)

        # Create a new parent clade and assign clades
        parent = Phylo.BaseTree.Clade(name=f"internal_{internal_counter}")
        internal_counter += 1
        parent.clades.append(clade1)
        parent.clades.append(clade2)

        clades.append(parent)

    tree = Phylo.BaseTree.Tree(root = clades[0])

    # Create a new root using normal_name and attach the built tree as its only child
    normal_root = Clade(name = normal_name)
    normal_root.clades.append(tree.root)
    tree.root = normal_root

    return tree


def _create_unique_random_start_trees(sample_labels, normal_name, n_chains, max_attempt_factor=20):
    """Create random start trees with unique topologies when possible."""
    trees = []
    topology_hashes = set()
    attempts = 0
    max_attempts = max_attempt_factor * max(1, n_chains)

    while len(trees) < n_chains and attempts < max_attempts:
        attempts += 1
        cur_tree = create_random_tree_topology(list(sample_labels), normal_name=normal_name)
        cur_hash = tree_hash.hash_tree_with_collapse(cur_tree)
        if cur_hash in topology_hashes:
            continue
        topology_hashes.add(cur_hash)
        trees.append(cur_tree)

    if len(trees) < n_chains:
        raise MEDICCError(
            "SPR mode: Could not generate unique random start trees for all chains. "
            f"Requested chains={n_chains}, unique_starts={len(trees)}, attempts={attempts}. "
            "Try reducing --n-cores or use --spr-start nj."
        )

    return trees


def update_branch_lengths(tree, fst, ancestor_fsa, normal_name='diploid'):
    """ Updates the branch lengths in the tree using the internal nodes supplied in the FSA dict 
    """
    if len(ancestor_fsa) == 2:
        child_clade = [x for x in tree.find_clades() if x.name is not None and x.name != normal_name][0]
        child_clade.branch_length = float(fstlib.score(
            fst, ancestor_fsa[normal_name], ancestor_fsa[child_clade.name]))

    if not isinstance(ancestor_fsa, dict):
        raise MEDICCError("input ancestor_fsa to function update_branch_lengths has to be either a dict"
                          "provided type is {}".format(type(ancestor_fsa)))

    def _distance_to_child(fst, ancestor_fsa, sample_1, sample_2):
        return float(fstlib.score(fst, ancestor_fsa[sample_1], ancestor_fsa[sample_2]))

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


def _spr_mode_single_chain(tree, samples_dict, fst, normal_name="diploid", prune_weight=0, MCMC_step=1000,
                           n_cores=None, sa_temp_start=10.0, sa_temp_end=0.01, sa_cooling="geometric",
                           seed=None, chain_id=1):
    """Perform tree search via SPR moves with simulated annealing.

    Uses incremental ancestor reconstruction: after each SPR move, only the
    up-pass (intersection) results for nodes on the prune/regraft-to-root paths
    are recomputed. The down-pass is always run in full.

    Simulated annealing schedule:
        - 'geometric': T(i) = T_start * alpha^i, where alpha is computed so that
          T(MCMC_step-1) = T_end. Standard for SPR-based phylogenetic search
          (used in e.g. GARLI).
        - 'linear': T(i) = T_start * (1 - i / MCMC_step). Simpler alternative.

    At each step, worse trees are accepted with probability exp(-delta / T(i)).
    High T → exploration (accept bad moves freely).
    Low T  → exploitation (mostly greedy, only accept improvements).

    Args:
        sa_temp_start: Starting temperature (default: 10.0).
        sa_temp_end: Final temperature (default: 0.01).
        sa_cooling: Cooling schedule, either 'geometric' or 'linear' (default: 'geometric').
    """
    from medicc.ancestors import reconstruct_ancestors_incremental

    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    # Compute cooling schedule parameters
    if sa_cooling == "geometric":
        # T(i) = T_start * alpha^i, solve for alpha such that T(N-1) = T_end
        if MCMC_step > 1:
            sa_alpha = (sa_temp_end / sa_temp_start) ** (1.0 / (MCMC_step - 1))
        else:
            sa_alpha = 1.0
        def temperature(step):
            return sa_temp_start * (sa_alpha ** step)
    elif sa_cooling == "linear":
        def temperature(step):
            return sa_temp_start * max(1e-10, 1.0 - step / MCMC_step)
    else:
        raise ValueError(f"Unknown cooling schedule '{sa_cooling}'. Use 'geometric' or 'linear'.")

    logger.info(
        f"SPR mode: chain {chain_id}, simulated annealing with {sa_cooling} cooling, "
        f"T_start={sa_temp_start}, T_end={sa_temp_end}, steps={MCMC_step}")

    sum_of_branch_length_l = []
    # Keep track of topologies visited
    topology_visited = set() # hash set used for proposing new topologies 
    full_tree_hash = set() # hash set used for keeping only true different topologies after branch length optimization, see comments in second part of tree_hash.py 

    # MCMC starting point — full reconstruction
    ancestors, uppass_cache = medicc.reconstruct_ancestors(tree=tree,
                                             samples_dict=samples_dict,
                                             fst=fst,
                                             normal_name=normal_name,
                                             prune_weight=prune_weight,
                                             spr_logger_disable=True,
                                             n_cores=n_cores)
    update_branch_lengths(tree, fst, ancestors, normal_name)
    sum_of_branch_length_l.append(medicc.tools.sum_of_branch_length(tree))

    topology_visited.add(tree_hash.tree_hash(tree_hash.strip_branch_lengths(tree)))
    full_tree_hash.add(tree_hash.hash_tree_with_collapse(tree)) 


    global_tree_l = [tree]
    global_ancestor_l = [ancestors]
    global_sum_of_branch_length = sum_of_branch_length_l[-1]

    tree_pre = tree
    sum_of_branch_length_pre = sum_of_branch_length_l[-1]
    ancestors_pre = ancestors
    uppass_cache_pre = uppass_cache
    for i in range(MCMC_step):
        T = temperature(i)

        # propose a new tree topology using SPR (with metadata for incremental reconstruction)
        spr_result = medicc.spr.spr_move_with_metadata(tree_pre)
        new_tree_topology = spr_result.tree

        # Check if the new topology has been visited before
        new_tree_topology_hash = tree_hash.tree_hash(tree_hash.strip_branch_lengths(new_tree_topology))
        if new_tree_topology_hash in topology_visited:
            logger.debug("SPR mode: step: {}, proposed new topology has been visited before, skipping".format(i))
            sum_of_branch_length_l.append(sum_of_branch_length_pre)
            continue
        else:
            topology_visited.add(new_tree_topology_hash)
            logger.debug("SPR mode: step: {}, proposed new topology is new".format(i))

        # Incremental ancestor reconstruction — only recompute dirty up-pass nodes
        # Pass the up-pass cache (pure intersection results), NOT the final ancestors
        new_ancestors, new_uppass_cache = reconstruct_ancestors_incremental(
            tree=new_tree_topology,
            old_uppass_cache=uppass_cache_pre,
            samples_dict=samples_dict,
            fst=fst,
            normal_name=normal_name,
            spr_result=spr_result,
            prune_weight=prune_weight)

        update_branch_lengths(new_tree_topology, fst, new_ancestors, normal_name)
        sum_of_branch_length_new = medicc.tools.sum_of_branch_length(new_tree_topology)
        logger.debug("SPR mode: step: {}, T: {:.4f}, proposed sum of branch length: {}".format(
            i, T, sum_of_branch_length_new))
            
        # Accept or reject the new tree (simulated annealing)
        if sum_of_branch_length_new <= sum_of_branch_length_pre:
            new_tree_topology_hash_full = tree_hash.hash_tree_with_collapse(new_tree_topology)
            if new_tree_topology_hash_full not in full_tree_hash:
                full_tree_hash.add(new_tree_topology_hash_full)
                logger.debug("SPR mode: step: {}, proposed new topology (after branch length optimization) is new and added to full_tree_hash".format(i))
                # Always accept improvements
                tree_pre = new_tree_topology
                sum_of_branch_length_pre = sum_of_branch_length_new
                ancestors_pre = new_ancestors
                uppass_cache_pre = new_uppass_cache
                sum_of_branch_length_l.append(sum_of_branch_length_new)
                logger.info("SPR mode: chain {}, step: {}, T: {:.4f}, accepted a better tree with sum of branch length {}".format(
                    chain_id, i, T, sum_of_branch_length_new))

                # Check if the new tree is better than the global best
                if sum_of_branch_length_new < global_sum_of_branch_length:
                    global_tree_l = [new_tree_topology]
                    global_ancestor_l = [new_ancestors]
                    global_sum_of_branch_length = sum_of_branch_length_new
                elif sum_of_branch_length_new == global_sum_of_branch_length:
                    global_tree_l.append(new_tree_topology)
                    global_ancestor_l.append(new_ancestors)
            else:
                logger.debug("SPR mode: step: {}, proposed new topology (after branch length optimization) is new but already in full_tree_hash, skipping".format(i))
                sum_of_branch_length_l.append(sum_of_branch_length_pre)
        else:
            # Accept worse tree with probability exp(-delta / T)
            delta = sum_of_branch_length_new - sum_of_branch_length_pre
            accept_prob = np.exp(-delta / T)
            if np.random.rand() < accept_prob:
                tree_pre = new_tree_topology
                sum_of_branch_length_pre = sum_of_branch_length_new
                ancestors_pre = new_ancestors
                uppass_cache_pre = new_uppass_cache
                sum_of_branch_length_l.append(sum_of_branch_length_new)
                logger.info("SPR mode: chain {}, step: {}, T: {:.4f}, accepted a worse tree (p={:.4f}) with sum of branch length {}".format(
                    chain_id, i, T, accept_prob, sum_of_branch_length_new))
            else:
                sum_of_branch_length_l.append(sum_of_branch_length_pre)
                logger.debug("SPR mode: step: {}, T: {:.4f}, rejected a worse tree (p={:.4f}) with sum of branch length {}".format(
                    i, T, accept_prob, sum_of_branch_length_new))

    return {
        "best_trees": global_tree_l,
        "best_ancestors": global_ancestor_l,
        "trace": sum_of_branch_length_l,
        "best_score": global_sum_of_branch_length,
    }


def spr_mode(tree, samples_dict, fst, normal_name="diploid", prune_weight=0, MCMC_step=1000,
             n_cores=None, sa_temp_start=10.0, sa_temp_end=0.01, sa_cooling="geometric",
             chain_start_trees=None):
    """Run SPR simulated annealing, using multi-chain mode automatically on multi-core runs."""
    if chain_start_trees is None:
        chain_start_trees = [copy.deepcopy(tree)]
    else:
        chain_start_trees = [copy.deepcopy(cur_tree) for cur_tree in chain_start_trees]

    if len(chain_start_trees) == 1:
        single_result = _spr_mode_single_chain(
            tree=chain_start_trees[0],
            samples_dict=samples_dict,
            fst=fst,
            normal_name=normal_name,
            prune_weight=prune_weight,
            MCMC_step=MCMC_step,
            n_cores=n_cores,
            sa_temp_start=sa_temp_start,
            sa_temp_end=sa_temp_end,
            sa_cooling=sa_cooling,
            chain_id=1,
        )
        return single_result["best_trees"], single_result["best_ancestors"], single_result["trace"]

    try:
        from joblib import Parallel, delayed
    except ImportError:
        raise ImportError("joblib must be installed for multi-chain SPR parallelization")

    n_chains = len(chain_start_trees)
    n_jobs = n_chains if n_cores is None else min(max(1, n_cores), n_chains)
    logger.info("SPR mode: Running %s parallel chains on %s cores", n_chains, n_jobs)

    seed_sequence = np.random.SeedSequence()
    chain_seeds = [int(seq.generate_state(1, dtype=np.uint32)[0]) for seq in seed_sequence.spawn(n_chains)]

    chain_results = Parallel(n_jobs=n_jobs)(
        delayed(_spr_mode_single_chain)(
            tree=chain_start_trees[chain_idx],
            samples_dict=samples_dict,
            fst=fst,
            normal_name=normal_name,
            prune_weight=prune_weight,
            MCMC_step=MCMC_step,
            n_cores=1,
            sa_temp_start=sa_temp_start,
            sa_temp_end=sa_temp_end,
            sa_cooling=sa_cooling,
            seed=chain_seeds[chain_idx],
            chain_id=chain_idx + 1,
        )
        for chain_idx in range(n_chains)
    )

    global_best_score = min(result["best_score"] for result in chain_results)
    global_tree_l = []
    global_ancestor_l = []
    seen_topology_hashes = set()

    for result in chain_results:
        if not np.isclose(result["best_score"], global_best_score):
            continue

        for cur_tree, cur_ancestors in zip(result["best_trees"], result["best_ancestors"]):
            cur_hash = tree_hash.hash_tree_with_collapse(cur_tree)
            if cur_hash in seen_topology_hashes:
                continue
            seen_topology_hashes.add(cur_hash)
            global_tree_l.append(cur_tree)
            global_ancestor_l.append(cur_ancestors)

    if len(global_tree_l) == 0:
        fallback_result = chain_results[0]
        global_tree_l = fallback_result["best_trees"]
        global_ancestor_l = fallback_result["best_ancestors"]

    chain_traces = [result["trace"] for result in chain_results]

    return global_tree_l, global_ancestor_l, chain_traces


def nni_mode(tree, samples_dict, fst, normal_name="diploid", prune_weight=0,
             nni_max_iter=1000, n_cores=None, nni_trace_dir=None):
    """Iterated steepest-ascent NNI hill-climbing with plateau traversal.

    At each sweep: enumerate all NNI neighbors of every tree in the current
    frontier, evaluate each via incremental ancestor reconstruction, and accept
    the globally best set:
      - If the best score strictly improves on current: replace the frontier
        with the single best neighbor (or all tied-best if multiple share it).
      - If the best score equals current: expand the frontier with all tied-best
        neighbors (plateau traversal — do not terminate yet).
      - If no neighbor improves or ties: terminate at an NNI-local optimum.

    Stop when a full sweep produces no improvement/tie, or nni_max_iter total
    neighbor trees have been evaluated.

    Args:
        nni_max_iter: Hard cap on the total number of neighbor trees evaluated
            across all sweeps (not the number of sweeps). Default 1000.
        nni_trace_dir: If not None, enables per-neighbor step recording. Used
            only as a boolean sentinel here — no files are written by nni_mode()
            itself; the caller (main_nni) is responsible for writing step_records
            to disk. Step numbers are pre-assigned in contiguous blocks before
            dispatch so parallel and serial paths produce identical numbering.

    Returns:
        dict with keys: best_trees (list of trees at optimum), best_ancestors
        (list of corresponding ancestor dicts), trace (list of scores, one
        entry per accepted iteration starting from initial), best_score (final
        score), step_records (list of (step, newick_str, score) tuples —
        empty when nni_trace_dir is None).
    """
    from medicc.ancestors import reconstruct_ancestors_incremental
    from medicc.nni import nni_neighbors, evaluate_nni_neighbors_parallel, _enumerate_moves
    from medicc.tree_hash import get_canonical_newick, get_topology_hash

    assert tree.root.name == normal_name, \
        "nni_mode expects a search-shape tree (root.name == normal_name)"

    # Initial full reconstruction
    ancestors, uppass_cache = medicc.reconstruct_ancestors(
        tree=tree, samples_dict=samples_dict, fst=fst,
        normal_name=normal_name, prune_weight=prune_weight,
        spr_logger_disable=True, n_cores=n_cores)
    update_branch_lengths(tree, fst, ancestors, normal_name)
    current_score = medicc.tools.sum_of_branch_length(tree)
    trace = [current_score]
    logger.info(f"NNI mode: initial score = {current_score}")

    # Frontier: list of (tree, ancestors, uppass_cache) tuples
    frontier = [(tree, ancestors, uppass_cache)]

    do_trace = nni_trace_dir is not None  # used as boolean sentinel; no files written here
    step_records = []  # list of (step, newick_str, score)
    step_offset = 0    # global across all sweeps — step numbers are unique for the entire run
    # step_offset doubles as the total-neighbors-evaluated counter for the nni_max_iter cap

    sweep = 0
    hit_cap = False
    frontier_hashes = frozenset(get_topology_hash(t) for t, _, _ in frontier)
    visited_hashes = set(frontier_hashes)  # all topology hashes ever in the frontier
    while True:
        best_score = None
        best_candidates = []  # list of (tree, ancestors, uppass_cache)
        n_evaluated = 0

        use_parallel = n_cores is not None and n_cores > 1
        for current_tree, current_ancestors, current_uppass_cache in frontier:
            if use_parallel:
                # Pre-assign step block before dispatch
                n_moves = len(_enumerate_moves(current_tree))
                step_start = step_offset
                step_offset += n_moves

                neighbor_results = evaluate_nni_neighbors_parallel(
                    tree=current_tree,
                    old_uppass_cache=current_uppass_cache,
                    samples_dict=samples_dict,
                    fst=fst,
                    normal_name=normal_name,
                    prune_weight=prune_weight,
                    n_cores=n_cores,
                    step_start=step_start,
                )
                for neighbor_tree, new_ancestors, new_uppass_cache, score, step in neighbor_results:
                    n_evaluated += 1
                    if do_trace:
                        step_records.append((step, get_canonical_newick(neighbor_tree), score))
                    if best_score is None or score < best_score:
                        best_score = score
                        best_candidates = [(neighbor_tree, new_ancestors, new_uppass_cache)]
                    elif score == best_score:
                        best_candidates.append((neighbor_tree, new_ancestors, new_uppass_cache))
            else:
                # Pre-assign step block before the serial loop
                n_moves = len(_enumerate_moves(current_tree))
                step_start = step_offset
                step_offset += n_moves

                for i, (neighbor_tree, synthetic_spr_result) in enumerate(nni_neighbors(current_tree)):
                    n_evaluated += 1
                    new_ancestors, new_uppass_cache = reconstruct_ancestors_incremental(
                        tree=neighbor_tree,
                        old_uppass_cache=current_uppass_cache,
                        samples_dict=samples_dict,
                        fst=fst,
                        normal_name=normal_name,
                        spr_result=synthetic_spr_result,
                        prune_weight=prune_weight)
                    update_branch_lengths(neighbor_tree, fst, new_ancestors, normal_name)
                    score = medicc.tools.sum_of_branch_length(neighbor_tree)
                    step = step_start + i
                    if do_trace:
                        step_records.append((step, get_canonical_newick(neighbor_tree), score))
                    if best_score is None or score < best_score:
                        best_score = score
                        best_candidates = [(neighbor_tree, new_ancestors, new_uppass_cache)]
                    elif score == best_score:
                        best_candidates.append((neighbor_tree, new_ancestors, new_uppass_cache))

        if not best_candidates:
            logger.info(f"NNI mode: sweep {sweep}: no neighbors enumerated, terminating")
            break

        def _dedup_candidates(candidates):
            seen = set()
            result = []
            for t, anc, cache in candidates:
                h = get_topology_hash(t)
                if h not in seen:
                    seen.add(h)
                    result.append((t, anc, cache))
            return result

        if best_score < current_score:
            previous_score = current_score
            current_score = best_score
            frontier = _dedup_candidates(best_candidates)
            frontier_hashes = frozenset(get_topology_hash(t) for t, _, _ in frontier)
            visited_hashes = set(frontier_hashes)  # reset visited on improvement
            trace.append(current_score)
            logger.info(
                f"NNI mode: sweep {sweep}: evaluated {n_evaluated} neighbors across "
                f"{len(frontier)} frontier tree(s), accepted improvement "
                f"{previous_score} -> {current_score} "
                f"({len(frontier)} unique tied-best neighbor(s), {len(best_candidates)} before dedup)")
        elif best_score == current_score:
            frontier = _dedup_candidates(best_candidates)
            new_frontier_hashes = frozenset(get_topology_hash(t) for t, _, _ in frontier)
            if new_frontier_hashes.issubset(visited_hashes):
                logger.info(
                    f"NNI mode: sweep {sweep}: evaluated {n_evaluated} neighbors across "
                    f"{len(frontier)} frontier tree(s), plateau exhausted — all "
                    f"{len(frontier)} neighbor(s) already visited at score {current_score}, terminating")
                break
            visited_hashes.update(new_frontier_hashes)
            frontier_hashes = new_frontier_hashes
            logger.info(
                f"NNI mode: sweep {sweep}: evaluated {n_evaluated} neighbors across "
                f"{len(frontier)} frontier tree(s), plateau — expanding frontier to "
                f"{len(frontier)} unique equally-good neighbor(s) at score {current_score} "
                f"({len(best_candidates)} before dedup)")
        else:
            logger.info(
                f"NNI mode: sweep {sweep}: evaluated {n_evaluated} neighbors across "
                f"{len(frontier)} frontier tree(s), no improvement "
                f"(best neighbor = {best_score}, current = {current_score}), "
                f"terminating at NNI-local optimum")
            break

        sweep += 1

        if step_offset >= nni_max_iter:
            hit_cap = True
            break

    if hit_cap:
        logger.warning(
            f"NNI mode: reached nni_max_iter={nni_max_iter} total neighbors evaluated "
            f"without converging — hard cap hit; returning current frontier of "
            f"{len(frontier)} tree(s) at score {current_score}")

    return {
        "best_trees": [t for t, _, _ in frontier],
        "best_ancestors": [a for _, a, _ in frontier],
        "trace": trace,
        "best_score": current_score,
        "step_records": step_records,
    }


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


def detect_wgd(input_df, sample, total_cn=False, wgd_x2=False, n_wgd=None):
    if n_wgd is not None and n_wgd > 2:
        raise NotImplementedError("MEDICC can only detect WGDs with n_wgd <= 2")

    if n_wgd is None:
        wgd_fst = io.read_fst(total_copy_numbers=total_cn, wgd_x2=wgd_x2, n_wgd=n_wgd)
        no_wgd_fst = io.read_fst(no_wgd=True)
    elif n_wgd == 1:
        wgd_fst = io.read_fst(total_copy_numbers=total_cn, wgd_x2=wgd_x2, n_wgd=2)
        no_wgd_fst = io.read_fst(total_copy_numbers=total_cn, wgd_x2=wgd_x2, n_wgd=1)
    elif n_wgd == 2:
        wgd_fst = io.read_fst(total_copy_numbers=total_cn, wgd_x2=wgd_x2, n_wgd=None)
        no_wgd_fst = io.read_fst(total_copy_numbers=total_cn, wgd_x2=wgd_x2, n_wgd=2)

    diploid_fsa = medicc.tools.create_diploid_fsa(no_wgd_fst)
    symbol_table = no_wgd_fst.input_symbols()
    fsa_dict, _ = medicc.create_standard_fsa_dict_from_data(input_df.loc[[sample]],
                                                         symbol_table, 'X')

    distance_wgd = float(fstlib.score(wgd_fst, diploid_fsa, fsa_dict[sample]))
    distance_no_wgd = float(fstlib.score(no_wgd_fst, diploid_fsa, fsa_dict[sample]))

    return distance_wgd < distance_no_wgd


class MEDICCError(Exception):
    pass
