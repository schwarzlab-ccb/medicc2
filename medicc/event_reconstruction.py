import logging
import os

import numpy as np
import pandas as pd

from medicc import io, native

# prepare logger 
logger = logging.getLogger(__name__)

EVENT_COLUMNS = ['sample_id', 'chrom', 'start', 'end', 'allele', 'type', 'cn_child']


def calculate_all_cn_events(tree, cur_df, alleles=['cn_a', 'cn_b'], normal_name='diploid',
                            model_config=None):
    """Create a DataFrame of minimum-cost copy-number events with the native model.

    For every branch one optimal ordered history (losses, gains and WGDs) is taken from the
    native engine and split into unit events on intervals of bins. Every event costs 1, so
    the number of events of a branch equals its branch length; the events replay the parent
    profile into the child profile in the order given by 'event_order'.

    Args:
        tree (Bio.Phylo.Tree): Final tree (rooted at the normal) with native branch lengths
        cur_df (pandas.DataFrame): Copy-numbers of the samples and internal nodes
        alleles (list): List of alleles
        normal_name (str, optional): Name of the normal sample. Defaults to 'diploid'.
        model_config (medicc.native.ModelConfig, optional): Native model. Defaults to the full model.

    Returns:
        pandas.DataFrame: Updated copy-number DataFrame (with is_gain, is_loss, is_wgd,
            is_normal and is_clonal)
        pandas.DataFrame: Events with the columns in EVENT_COLUMNS plus 'event_order' (order of
            the event on its branch) and 'stage' (index of the history stage)

    Raises:
        MEDICCEventReconstructionError: If a history does not replay or its number of events
            differs from the branch length.
    """
    from medicc.core import create_cn_str_dict

    cur_df = cur_df.copy()
    cur_df[alleles] = cur_df[alleles].astype(int)
    cur_df[['is_gain', 'is_loss', 'is_wgd']] = False
    profiles = create_cn_str_dict(cur_df[alleles])
    # (allele, coordinates of one chromosome) in the order of the X-separated blocks
    normal_index = cur_df.loc[normal_name].index
    chromosomes = [list(x.index) for _, x in
                   pd.Series(0, index=normal_index).groupby(level='chrom', observed=False)]
    blocks = [(allele, coords) for allele in alleles for coords in chromosomes]
    genome_start = int(normal_index.get_level_values('start').min())
    genome_end = int(normal_index.get_level_values('end').max())

    events = []
    with native.Model(model_config) as model:
        for clade in tree.find_clades():
            parent_name = clade.name if clade.name is not None else normal_name
            for child in clade.clades:
                if child.name == normal_name:
                    continue
                score, stages = model.history(profiles[parent_name], profiles[child.name])
                if score != child.branch_length:
                    raise MEDICCEventReconstructionError(
                        f"Event history score {score} differs from the branch length "
                        f"{child.branch_length} of {child.name}")
                branch_events = []
                for stage_nr, stage in enumerate(stages):
                    if stage['type'] == 'wgd':
                        doubled = ''.join(c if c == 'X' else str(min(8, 2 * int(c)))
                                          for c in stage['before'])
                        if doubled != stage['after']:
                            raise MEDICCEventReconstructionError("Invalid WGD in native history")
                        branch_events.append([child.name, 'chr0', genome_start, genome_end,
                                              'both', 'wgd', 0, len(branch_events) + 1, stage_nr])
                        cur_df.loc[child.name, 'is_wgd'] = True
                        continue
                    for (allele, coords), before, after in zip(
                            blocks, stage['before'].split('X'), stage['after'].split('X')):
                        branch_events.extend(_native_stage_events(
                            cur_df, child.name, allele, coords, before, after, stage['type'],
                            stage_nr, len(branch_events)))
                if len(branch_events) != score:
                    raise MEDICCEventReconstructionError(
                        f"Number of events ({len(branch_events)}) differs from the branch length "
                        f"({score}) of {child.name}")
                events.extend(branch_events)

    events_df = pd.DataFrame(events, columns=EVENT_COLUMNS + ['event_order', 'stage'])
    events_df = events_df.set_index(['sample_id', 'chrom', 'start', 'end'])

    flags = cur_df[['is_gain', 'is_loss', 'is_wgd']].any(axis=1)
    is_normal = ~flags.groupby(level=['chrom', 'start', 'end'], observed=True).any()
    mrca = [x for x in tree.root.clades if x.name != normal_name][0].name
    is_clonal = ~flags.loc[flags.index.get_level_values('sample_id') != mrca].groupby(
        level=['chrom', 'start', 'end'], observed=True).any()
    positions = cur_df.index.droplevel('sample_id')
    cur_df['is_normal'] = is_normal.reindex(positions).to_numpy()
    cur_df['is_clonal'] = is_clonal.reindex(positions).to_numpy()

    return cur_df, events_df


def _native_stage_events(cur_df, sample_id, allele, coords, before, after, event_type, stage_nr,
                         nr_previous):
    """Split one loss or gain stage of one chromosome into unit interval events and mark the
    affected bins in cur_df. Returns the event rows (EVENT_COLUMNS, event_order, stage)."""
    before = np.array(list(before), dtype=int)
    after = np.array(list(after), dtype=int)
    current = before.copy()
    rows = []
    for first, end in _unit_intervals(before, after, event_type):
        old = current[first:end].copy()
        if event_type == 'loss':
            current[first:end] = np.maximum(0, old - 1)
        else:
            current[first:end] = np.where(old == 0, 0, old + 1)
        for i in range(first, end):
            if old[i - first] != current[i]:
                cur_df.loc[(sample_id, *coords[i]), 'is_' + event_type] = True
        rows.append([sample_id, coords[first][0], coords[first][1], coords[end - 1][2], allele,
                     event_type, int(current[first]), nr_previous + len(rows) + 1, stage_nr])
    if not np.array_equal(current, after):
        raise MEDICCEventReconstructionError("Interval events do not replay the native history")
    return rows


def _unit_intervals(before, after, event_type):
    """Minimal set of unit loss or gain intervals (start, end) that turns before into after.
    Bins that are zero after the stage can be bridged by an interval at no extra cost."""
    delta = before - after if event_type == 'loss' else after - before
    if (delta < 0).any():
        raise MEDICCEventReconstructionError("History stage has the wrong direction")
    bridge = after == 0
    for height in range(1, int(delta.max(initial=0)) + 1):
        required = delta >= height
        first = last = None
        for i in range(len(delta) + 1):
            if i < len(delta) and required[i]:
                if first is None:
                    first = i
                last = i
            elif i == len(delta) or not bridge[i]:
                if first is not None:
                    yield first, last + 1
                    first = last = None


def overlap_events(events_df=None, output_df=None, tree=None, overlap_threshold=0.9,
                   chromosome_bed='default', regions_bed='default',
                   replace_loh_with_loss=True, alleles=['cn_a', 'cn_b'],
                   replace_both_arms_with_chrom=True, normal_name='diploid'):
    """Overlap copy-number events with regions of interest

    Args:
        events_df (pandas.DataFrame, optional): All copy-number events. Defaults to None.
        output_df (pandas.DataFrame, optional): DataFrame containing all copy-numbers. Defaults to None.
        tree (Bio.Phylo.Tree, optional): Phylogenetic tree. Defaults to None.
        overlap_threshold (float, optional): Threshold above which an overlap is considered. Defaults to 0.9.
        chromosome_bed (str, optional): Name of BED file containing chromosome arm information. Defaults to 'none'.
        regions_bed (str, optional): Name of BED file containing regions of interest. Defaults to 'none'.
        replace_loh_with_loss (bool, optional): If True, loh is considered like a normal loss. Defaults to True.
        alleles (list, optional): List of alleles. Defaults to ['cn_a', 'cn_b'].
        replace_both_arms_with_chrom (bool, optional): If True, an event in the p- and q-arm of a chromosome will be displayed as a single event. Defaults to True.
        normal_name (str, optional): Name of normal sample. Defaults to 'diploid'.

    Returns:
        pandas.DataFrame: DataFrame with events concerning the regions of interest
    """
    try:
        import pyranges as pr
    except ImportError:
        raise ImportError("You have to install pyranges to overlap events with regions of interest") from None       

    if chromosome_bed == 'none':
        chromosome_bed = None
    elif chromosome_bed == 'default':
        chromosome_bed = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                      "objects", "hg38_chromosome_arms.bed")
    if regions_bed == 'none':
        regions_bed = None
    elif regions_bed == 'default':
        regions_bed = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   "objects", "Davoli_2013_TSG_OG_genes.bed")

    all_events = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'name', 'NumberOverlaps',
                                       'FractionOverlaps', 'event', 'branch']).set_index(['Chromosome', 'Start', 'End'])

    if events_df is None:
        if output_df is None or tree is None:
            raise MEDICCEventReconstructionError("Either events_df or df and tree has to be specified")
        _, events_df = calculate_all_cn_events(tree, output_df, alleles=alleles, normal_name=normal_name)
    if replace_loh_with_loss:
        events_df.loc[events_df['type'] == 'loh', 'type'] = 'loss'

    # Read chromosome regions and other regions
    if chromosome_bed is None and regions_bed is None:
        raise MEDICCEventReconstructionError("Either chromosome_bed or regions_bed has to be specified")

    chr_arm_regions = None
    if chromosome_bed is not None:
        logger.debug(f'Overlap with chromosomes bed file {chromosome_bed}')
        chr_arm_regions = io.read_bed_file(chromosome_bed)
        whole_chromosome = chr_arm_regions.groupby('Chromosome').min()
        whole_chromosome['End'] = chr_arm_regions.groupby('Chromosome')['End'].max()
        whole_chromosome['name'] = whole_chromosome.index
        chr_arm_regions = pd.concat([chr_arm_regions, whole_chromosome.reset_index()]).sort_values('Chromosome')
        chr_arm_regions = pr.PyRanges(chr_arm_regions)

    regions = None
    if regions_bed is not None:
        logger.debug(f'Overlap with regions bed file {regions_bed}')
        regions = []
        if isinstance(regions_bed, list) or isinstance(regions_bed, tuple):
            for f in regions_bed:
                regions.append(pr.PyRanges(io.read_bed_file(f)))
        else:
            regions.append(pr.PyRanges(io.read_bed_file(regions_bed)))

    # Add WGD: one row per WGD event
    wgd_branches = events_df.loc[events_df['type'] == 'wgd'].index.get_level_values('sample_id')
    wgd_index = pd.MultiIndex.from_tuples([('all', 0, 0)] * len(wgd_branches),
                                          names=['Chromosome', 'Start', 'End'])
    wgd_events = pd.DataFrame({'name': 'WGD', 'NumberOverlaps': 1, 'FractionOverlaps': 1.,
                               'event': 'WGD', 'branch': list(wgd_branches)},
                              index=wgd_index, columns=all_events.columns)

    overlaps = [all_events, wgd_events]
    for cur_branch in events_df.index.get_level_values('sample_id').unique():
        for allele in alleles:
            cur_events_df = events_df.loc[cur_branch]
            cur_events_df = cur_events_df.loc[cur_events_df['allele']==allele]
            for event_type in ['gain', 'loss'] if replace_loh_with_loss else ['gain', 'loh', 'loss']:
                cur_events_ranges = pr.PyRanges(cur_events_df.loc[cur_events_df['type'] == event_type].reset_index(
                ).rename({'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}, axis=1))

                # Calculate chromosomal events
                if chr_arm_regions is not None:
                    chr_events = overlap_regions(
                        chr_arm_regions, cur_events_ranges, event_type, cur_branch, overlap_threshold)
                    # remove arms if the whole chromosome is in there
                    if replace_both_arms_with_chrom and len(chr_events) > 0:
                        chr_events = chr_events[~chr_events['name'].isin(np.concatenate(
                            [[name + 'p', name + 'q'] if ('q' not in name and 'p' not in name) else [] for name in chr_events['name']]))]
                    overlaps.append(chr_events)

                # Calculate other events
                if regions is not None:
                    for region in regions:
                        overlaps.append(overlap_regions(
                            region, cur_events_ranges, event_type, cur_branch, overlap_threshold))

    overlaps = [x for x in overlaps if len(x) > 0]
    if len(overlaps) > 0:
        all_events = pd.concat(overlaps)

    all_events['final_name'] = all_events['name'].apply(lambda x: x.split(
        'chr')[-1]) + all_events['event'].apply(lambda x: ' +' if x == 'gain' else (' -' if x == 'loss' else (' 0' if x == 'loh' else '')))

    all_events = (all_events
                  .reset_index()
                  .sort_values(['branch', 'Chromosome', 'Start', 'End'], ascending=True)
                  .set_index(['branch'])
                  .drop('NumberOverlaps', axis=1)
                  .rename({'Start': 'Event-Start', 'End': 'Event-End'}, axis=1)
                  )

    return all_events


def overlap_regions(region, cur_events_ranges, event, branch, overlap_threshold):

    cur_events_overlaps = region.coverage(cur_events_ranges).as_df()
    cur_events_overlaps = cur_events_overlaps.loc[cur_events_overlaps['FractionOverlaps']
                                                > overlap_threshold]
    cur_events_overlaps = cur_events_overlaps.set_index(['Chromosome', 'Start', 'End'])
    cur_events_overlaps['event'] = event
    cur_events_overlaps['branch'] = branch

    return cur_events_overlaps


class MEDICCEventReconstructionError(Exception):
    pass
