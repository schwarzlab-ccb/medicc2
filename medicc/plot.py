import logging

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

COL_ALLELE_A = mpl.colors.to_rgba('orange')
COL_ALLELE_B = mpl.colors.to_rgba('teal')
COL_CLONAL = mpl.colors.to_rgba('lightgrey')
COL_NORMAL = mpl.colors.to_rgba('dimgray')
COL_GAIN = mpl.colors.to_rgba('red')
COL_WGD = mpl.colors.to_rgba('green')
COL_LOSS = mpl.colors.to_rgba('blue')
COL_CHR_LABEL = mpl.colors.to_rgba('grey')
COL_VLINES = '#1f77b4'
COL_MARKER_INTERNAL = COL_VLINES
COL_MARKER_TERMINAL = 'black'
COL_MARKER_NORMAL = 'green'
COL_SUMMARY_LABEL = 'grey'
COL_BACKGROUND = 'white'
COL_BACKGROUND_HATCH = 'lightgray'
COL_PATCH_BACKGROUND = 'white'
LINEWIDTH_COPY_NUMBERS = 2
LINEWIDTH_CHR_BOUNDARY = 1
LINEWIDTH_SEGMENT_BOUNDARY = 0.5
ALPHA_PATCHES = 0.15
ALPHA_PATCHES_WGD = 0.3
ALPHA_CLONAL = 0.3
BACKGROUND_HATCH_MARKER = '/////'
TREE_MARKER_SIZE = 40
YLABEL_FONT_SIZE = 8
YLABEL_TICK_SIZE = 6
XLABEL_FONT_SIZE = 10
XLABEL_TICK_SIZE = 8
CHR_LABEL_SIZE = 8
SMALL_SEGMENTS_LIMIT = 1e7


def plot_cn_profiles(
        input_df,
        input_tree=None,
        title=None,
        normal_name='diploid',
        mincn='auto',
        maxcn='auto',
        allele_columns=['cn_a', 'cn_b'],
        plot_summary=True,
        plot_subclonal_summary=True,
        plot_clonal_summary=False,
        hide_normal_chromosomes=False,
        ignore_segment_lengths=False,
        tree_width_scale=1,
        track_width_scale=1,
        height_scale=1,
        horizontal_margin_adjustment=0,
        close_gaps=False,
        show_small_segments=False,
        show_branch_support=False,
        hide_internal_nodes=False,
        chr_label_func=None,
        show_events_in_tree=False,
        show_events_in_cn=True,
        show_branch_lengths=True,
        detailed_xticks=False,
        clonal_transparant=False,
        label_func=None):
    
    df = input_df.copy()
    df[allele_columns] = df[allele_columns].astype(int)
    if len(allele_columns) > 2:
        logger.warning("More than two allels were provided ({})\n"
                    "Copy number tracks can only be plotted for 1 or 2 alleles".format(allele_columns))
    if len(np.setdiff1d(allele_columns, df.columns)):
        logger.warning("Some provided allele_columns are not in the dataframe"
                    "These are: {}".format(np.setdiff1d(allele_columns, df.columns)))

    if np.setdiff1d(['is_clonal', 'is_normal', 'is_gain', 'is_loss', 'is_wgd'], df.columns).size > 0:
        if input_tree is None:
            df[['is_normal', 'is_clonal', 'is_gain', 'is_loss', 'is_wgd']] = False
        else:
            df[['is_gain', 'is_loss', 'is_wgd']] = False

            cn_change = compute_cn_change(df=df[allele_columns], tree=input_tree, normal_name=normal_name)
            df.loc[(cn_change < 0).any(axis=1), 'is_loss'] = True
            df.loc[(cn_change > 0).any(axis=1), 'is_gain'] = True

            is_normal = ~df.unstack('sample_id')[['is_loss', 'is_gain', 'is_wgd']].any(axis=1)
            is_normal.name = 'is_normal'
            mrca = [x for x in input_tree.root.clades if x.name != normal_name][0].name
            is_clonal = ~df.loc[df.index.get_level_values('sample_id')!=mrca].unstack('sample_id')[['is_loss', 'is_gain', 'is_wgd']].any(axis=1)
            is_clonal.name = 'is_clonal'

            df = (df
                  .drop(['is_normal', 'is_clonal'], axis=1, errors='ignore')
                  .join(is_normal, on=['chrom', 'start', 'end'], how='inner')
                  .join(is_clonal, on=['chrom', 'start', 'end'], how='inner')
                  .sort_index()
                  )

    if hide_normal_chromosomes:
        df = df.join(df.groupby('chrom')['is_normal'].all().to_frame('hide'))
        df = df.query("~hide").drop('hide', axis=1)

    if mincn=='auto':
        mincn = df[allele_columns].min().min()
    else:
        mincn = int(mincn)        
    if maxcn=='auto':
        maxcn = df[allele_columns].max().max()
    else:
        maxcn = int(maxcn)

    samples = df.index.get_level_values('sample_id').unique()
    if hide_internal_nodes:
        if input_tree is not None:
            samples = [x for x in df.index.get_level_values('sample_id').unique() if list(
                input_tree.find_clades(x))[0].is_terminal()]
        else:
            samples = [x for x in df.index.get_level_values('sample_id').unique() if 'internal_' not in x]
    nsamp = len(samples)
    nsegs = df.loc[samples[0], :].groupby('chrom').size()

    # == 2 because of diploid
    if input_tree is None or normal_name is None or nsamp == 2:
        plot_summary = False
        plot_subclonal_summary = False
        plot_clonal_summary = False

    if nsamp > 20:
        logger.warning('More than 20 samples were provided. Creating the copy number tracks will take '
                    'a long time to process and might crash. Best to use plot_cn_heatmap instead')

    df.reset_index(['start','end'], inplace=True)

    segment_lengths = df['end'] - df['start']
    df['small_segment'] = (segment_lengths) < SMALL_SEGMENTS_LIMIT
    if ignore_segment_lengths:
        segment_lengths.loc[:] = 1

    if close_gaps or ignore_segment_lengths:
        cur_df = df.loc[samples[0], :].reset_index()[['chrom', 'start', 'end']]

        cur_df['start_pos'] = np.cumsum(np.append([0], segment_lengths.loc[samples[0]]))[:-1]
        cur_df['end_pos'] = cur_df['start_pos'] + segment_lengths.loc[samples[0]].values

        cur_df.set_index(['chrom', 'start'], inplace=True)
        cur_df.drop('end', axis=1, inplace=True)

        df = df.join(cur_df, on=['chrom', 'start'])

    else:
        chrom_offset = df.loc[samples[0],:].reset_index().groupby('chrom', sort=False).max()['end']
        chrom_offset.dropna(inplace=True)
        chrom_offset.loc[:] = np.append(0, chrom_offset.cumsum().values[:-1])
        chrom_offset.name='offset'
        df = df.join(chrom_offset, on='chrom')
        df['start_pos'] = df['start'] + df['offset']
        df['end_pos'] = df['end'] + df['offset']

    df = df.reset_index().set_index(['sample_id', 'chrom', 'start', 'end'])

    if plot_summary or plot_subclonal_summary:
        agg_events = compute_cn_change(df=df[allele_columns], tree=input_tree, normal_name=normal_name)
        agg_events = agg_events.groupby(["chrom", "start", "end"], observed=True).sum()

        agg_events[['start_pos', 'end_pos', 'small_segment']
                       ] = df.loc[samples[0], ['start_pos', 'end_pos', 'small_segment']]

        mrca = [x for x in input_tree.root.clades if x.name != normal_name][0].name
        mrca_df = df.loc[mrca].copy()
        mrca_df.loc[:, allele_columns] = mrca_df.loc[:, allele_columns] - df.loc[normal_name, allele_columns]

    ## determine clade colors
    clade_colors = {}
    for sample in samples:
        ## determine if sample is terminal
        is_terminal = True
        if input_tree is not None:
            matches = list(input_tree.find_clades(sample))
            if len(matches)>0:
                clade = matches[0]
                is_terminal = clade.is_terminal()
        ## determine if sample is normal
        clade_colors[sample] = COL_MARKER_TERMINAL
        if not is_terminal:
            clade_colors[sample] = COL_MARKER_INTERNAL
        if sample == normal_name:
            clade_colors[sample] = COL_MARKER_NORMAL

    # define plot width and height and create figure
    track_width = nsegs.sum() * 0.2 * track_width_scale
    if input_tree is not None:
        tree_width = 2.5 * tree_width_scale ## in figsize
    else:
        tree_width = 0

    nrows = nsamp + int(plot_summary) + int(plot_subclonal_summary) + int(plot_clonal_summary)
    plotheight = 4 * 0.2 * nrows * height_scale
    plotwidth = tree_width + track_width
    tree_width_ratio = tree_width / plotwidth
    fig = plt.figure(figsize=(min(250, plotwidth), min(250, plotheight)), constrained_layout=True)
    if input_tree is None:
        gs = fig.add_gridspec(nrows, 1)
        cn_axes = [fig.add_subplot(gs[i]) for i in range(0, nrows)]
        y_posns = {sample: i for i, sample in enumerate(samples)} # as they appear
    else:
        gs = fig.add_gridspec(nrows, 2, width_ratios=[tree_width_ratio, 1-tree_width_ratio])
        tree_ax = fig.add_subplot(gs[0:nsamp, 0])
        cn_axes = [fig.add_subplot(gs[i]) for i in range(1,(2*(nrows))+1,2)]
        y_posns = _get_y_positions(
            input_tree, adjust=not hide_internal_nodes, normal_name=normal_name)
        y_posns = {clade.name: y_pos for clade, y_pos in y_posns.items(
            ) if clade.name is not None and clade.name != 'root'}
        plot_tree(input_tree, 
                  ax=tree_ax,
                  title=title,
                  normal_name=normal_name,
                  label_func=lambda x: '',
                  label_colors=clade_colors,
                  show_branch_support=show_branch_support,
                  show_events=show_events_in_tree,
                  show_branch_lengths=show_branch_lengths,
                  hide_internal_nodes=hide_internal_nodes)
    
    # Adjust the margin between the tree and cn tracks. Default is 0
    fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=horizontal_margin_adjustment)
    ## iterate over samples and plot the track
    for sample, group in df.groupby('sample_id'):
        if sample not in samples:
            continue
        index_to_plot = y_posns[sample] - 1
        plot_axis_labels = (index_to_plot == (nsamp-1))
        _plot_cn_profile(ax=cn_axes[index_to_plot],
                         label=label_func(sample) if label_func is not None else sample,
                         data=group,
                         mincn=mincn-1,
                         maxcn=maxcn+1,
                         alleles=allele_columns,
                         type='sample',
                         clonal_transparant=clonal_transparant,
                         chr_label_func=chr_label_func,
                         plot_xaxis_labels=plot_axis_labels if not (
                            plot_summary + plot_subclonal_summary + plot_clonal_summary) else False,
                         plot_yaxis_labels=True,
                         yaxis_label_color=clade_colors[sample],
                         show_small_segments=show_small_segments,
                         detailed_xticks=detailed_xticks,
                         show_events_in_cn=show_events_in_cn)

    if plot_clonal_summary:
        cur = mrca_df.copy()
        cur[['is_clonal', 'is_normal', 'is_gain', 'is_loss', 'is_wgd']] = False
        _plot_cn_profile(ax=cn_axes[nsamp],
                         label='clonal\nchanges',
                         data=cur,
                         mincn=mrca_df[allele_columns].min().min()-1,
                         maxcn=mrca_df[allele_columns].max().max()+1,
                         alleles=allele_columns,
                         detailed_xticks=detailed_xticks,
                         chr_label_func=chr_label_func,
                         type='summary',
                         show_small_segments=show_small_segments)
        cn_axes[nsamp].get_xaxis().set_visible(not (plot_summary or plot_subclonal_summary))

    if plot_summary:
        cur = agg_events.copy()
        cur[['is_clonal', 'is_normal', 'is_gain', 'is_loss', 'is_wgd']] = False
        _plot_cn_profile(ax=cn_axes[-1],
                         label='all\nchanges',
                         data=cur,
                         mincn=agg_events[allele_columns].min().min()-1,
                         maxcn=agg_events[allele_columns].max().max()+1,
                         alleles=allele_columns,
                         detailed_xticks=detailed_xticks,
                         chr_label_func=chr_label_func,
                         type='summary',
                         show_small_segments=show_small_segments)

    if plot_subclonal_summary:
        agg_events.loc[:, allele_columns] = agg_events[allele_columns] - mrca_df[allele_columns]
        cur = agg_events.copy()
        cur[['is_clonal', 'is_normal', 'is_gain', 'is_loss', 'is_wgd']] = False

        _plot_cn_profile(ax=cn_axes[-1 - int(plot_summary)],
                         label='subclonal\nchanges',
                         data=cur,
                         type='summary',
                         mincn=agg_events[allele_columns].min().min()-1,
                         maxcn=agg_events[allele_columns].max().max()+1,
                         alleles=allele_columns,
                         detailed_xticks=detailed_xticks,
                         chr_label_func=chr_label_func,
                         show_small_segments=show_small_segments)
        cn_axes[-1 - int(plot_summary)].get_xaxis().set_visible(not plot_summary)

    cn_axes[-1].set_xlabel('genome position', fontsize=XLABEL_FONT_SIZE)
    cn_axes[-1].xaxis.set_tick_params(labelsize=XLABEL_TICK_SIZE)
    cn_axes[-1].xaxis.offsetText.set_fontsize(XLABEL_TICK_SIZE)

    return fig


def _plot_cn_profile(ax, label, data, mincn, maxcn, alleles, type='sample',
                     plot_xaxis_labels=True, plot_yaxis_labels=True, chr_label_func=None,
                     clonal_transparant=True, yaxis_label_color='black', show_small_segments=False,
                     show_events_in_cn=True, detailed_xticks=False):
    ## collect line segments and background patches
    event_patches = []
    bkg_patches = []
    lines_a = []
    lines_b = []
    circles_a = []
    circles_b = []

    two_alleles = len(alleles) == 2

    alpha = []
    for idx, r in data.iterrows():
        lines_a.append([(r['start_pos'], r[alleles[0]]), 
                        (r['end_pos'], r[alleles[0]])])
        if two_alleles:
            lines_b.append([(r['start_pos'], r[alleles[1]]), (r['end_pos'], r[alleles[1]])])
        rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'],
                                     maxcn-mincn, edgecolor=None, facecolor=COL_PATCH_BACKGROUND, alpha=1)
        alpha.append(ALPHA_CLONAL if (clonal_transparant and r['is_clonal'] and type=='sample') else 1.0)
        bkg_patches.append(rect)

        if show_small_segments and r['small_segment']:
            circles_a.append((r['start_pos'] + 0.5*(r['end_pos'] - r['start_pos']), r[alleles[0]]))
            if two_alleles:
                circles_b.append(
                    (r['start_pos'] + 0.5*(r['end_pos'] - r['start_pos']), r[alleles[1]]))

        if show_events_in_cn and type=='sample':
            if r['is_normal']:
                rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'],
                                            maxcn-mincn, edgecolor=None, facecolor=COL_NORMAL, alpha=ALPHA_PATCHES)
                event_patches.append(rect)
            if r['is_wgd']:
                rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'],
                                             maxcn-mincn, edgecolor=None, facecolor=COL_WGD, alpha=ALPHA_PATCHES_WGD)
                event_patches.append(rect)
            if r['is_gain']:
                rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'],
                                            maxcn-mincn, edgecolor=None, facecolor=COL_GAIN, alpha=ALPHA_PATCHES)
                event_patches.append(rect)
            if r['is_loss']:
                rect = mpl.patches.Rectangle((r['start_pos'], mincn), r['end_pos']-r['start_pos'],
                                            maxcn-mincn, edgecolor=None, facecolor=COL_LOSS, alpha=ALPHA_PATCHES)
                event_patches.append(rect)

    events = mpl.collections.PatchCollection(event_patches, match_original=True, zorder=2)
    ax.add_collection(events)
    backgrounds = mpl.collections.PatchCollection(bkg_patches, match_original=True, zorder=1)
    ax.add_collection(backgrounds)
    plot_bkg = mpl.patches.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=COL_BACKGROUND,
                                     edgecolor=COL_BACKGROUND_HATCH, zorder=0, hatch=BACKGROUND_HATCH_MARKER)
    ax.add_patch(plot_bkg)
    colors_a = np.array([COL_ALLELE_A] * len(lines_a))
    # clonal mutations
    colors_a[:, 3] = np.array(alpha)
    if two_alleles:
        colors_b = np.array([COL_ALLELE_B] * len(lines_b))
        colors_b[:, 3] = np.array(alpha)
        # a and b are overlapping
        colors_a[data[alleles[0]] == data[alleles[1]], 3] = 0.5
        colors_b[data[alleles[0]] == data[alleles[1]], 3] = 0.5
        colors = np.vstack([colors_a, colors_b])
    else:
        colors = colors_a
    lc = mpl.collections.LineCollection(
        lines_a + lines_b, colors=colors, linewidth=LINEWIDTH_COPY_NUMBERS)
    ax.add_collection(lc)

    if type == 'summary':
        ax.axhline(0, color='grey', linestyle='--', zorder=1,
                   linewidth=LINEWIDTH_SEGMENT_BOUNDARY, alpha=0.75,)

    if len(circles_a) > 0:
        if two_alleles:
            a_b_overlap = (data.loc[data['small_segment']][alleles[0]] ==
                           data.loc[data['small_segment']][alleles[1]]).values
        else:
            a_b_overlap = circles_a * [False]
        ax.plot(np.array(circles_a)[~a_b_overlap, 0], np.array(circles_a)[~a_b_overlap, 1],
                'o', ms=3, color=COL_ALLELE_A, alpha=1., zorder=6)
        ax.plot(np.array(circles_a)[a_b_overlap, 0], np.array(circles_a)[a_b_overlap, 1],
                'o', ms=3, color=COL_ALLELE_A, alpha=0.5, zorder=6)
        if two_alleles:
            ax.plot(np.array(circles_b)[~a_b_overlap, 0], np.array(circles_b)[~a_b_overlap, 1],
                    'o', ms=3, color=COL_ALLELE_B, alpha=1., zorder=6)
            ax.plot(np.array(circles_b)[a_b_overlap, 0], np.array(circles_b)[a_b_overlap, 1],
                    'o', ms=3, color=COL_ALLELE_B, alpha=0.5, zorder=6)

    ax.autoscale()

    ## draw segment boundaries
    seg_bound_first = data['start_pos'].values[0]
    seg_bounds = data['end_pos'].values
    ax.vlines(np.append(seg_bound_first, seg_bounds), ymin=mincn, ymax=maxcn, ls='--', alpha=0.25,
              color=COL_VLINES, linewidth=LINEWIDTH_SEGMENT_BOUNDARY)

    ## draw chromosome boundaries
    chr_ends = data.reset_index().groupby('chrom').max()['end_pos']
    chr_ends.dropna(inplace=True)
    linex = chr_ends.values[:-1]  # don't plot last
    ax.vlines(linex, ymin=mincn, ymax=maxcn, color=COL_VLINES, linewidth=LINEWIDTH_CHR_BOUNDARY)

    ## draw chromosome labels
    chr_label_pos = chr_ends
    chr_label_pos.loc[:] = np.roll(chr_label_pos.values, 1)
    chr_label_pos.iloc[0] = seg_bound_first
    for chrom, pos in chr_label_pos.items():
        if chr_label_func is not None:
            chromtxt = chr_label_func(chrom)
        else:
            chromtxt = chrom

        ax.text(pos, maxcn-0.35, chromtxt, ha='left', va='top', color=COL_CHR_LABEL,
                fontweight='medium', fontsize=CHR_LABEL_SIZE)

    ## draw sample labels
    if plot_yaxis_labels:
        ax.set_ylabel(label, fontsize=YLABEL_FONT_SIZE, rotation=0, ha='right', va='center')
        ax.yaxis.label.set_color(yaxis_label_color)

    ## axis modifications
    if detailed_xticks:
        ax.set_xticks(np.arange(0, data['end_pos'].max(), 1e7))
    ax.get_xaxis().set_visible(plot_xaxis_labels)

    #ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune='both', nbins=(maxcn-mincn)/2 + 1))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune='both', nbins=5))
    ax.yaxis.set_tick_params(labelsize=YLABEL_TICK_SIZE)
    ax.set_ylim(mincn, maxcn)
    ax.set_xlim(0, data['end_pos'].max())

def _get_x_positions(tree):
    """Create a mapping of each clade to its horizontal position.
    Dict of {clade: x-coord}
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths

def _get_y_positions(tree, adjust=False, normal_name='diploid'):
    """Create a mapping of each clade to its vertical position.
    Dict of {clade: y-coord}.
Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    heights = {tip: maxheight -1 -i for i,
            tip in enumerate(reversed([x for x in tree.get_terminals() if x.name != normal_name]))}
    if normal_name is not None:
        heights.update({list(tree.find_clades(normal_name))[0]: maxheight})

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
        
    if adjust:
        pos = pd.DataFrame([(clade, val) for clade, val in heights.items()], columns=['clade','pos']).sort_values('pos')
        pos['newpos'] = 0
        count = 0
        for i in pos.index:
            if pos.loc[i,'clade'].name is not None and pos.loc[i,'clade'].name != 'root':
                count = count+1
            pos.loc[i, 'newpos'] = count

        pos.set_index('clade', inplace=True)
        heights = pos.to_dict()['newpos']

    return heights


def plot_tree(input_tree,
              label_func=None,
              title='',
              ax=None,
              output_name=None,
              normal_name='diploid',
              width_scale=1,
              height_scale=1,
              show_branch_lengths=True,
              show_branch_support=False,
              show_events=False,
              branch_labels=None,
              label_colors=None,
              hide_internal_nodes=False,
              marker_size=None,
              line_width=None,
              **kwargs):
    """Plot the given tree using matplotlib (or pylab).
    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.
    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).
    Example using the pyplot options 'axhspan' and 'axvline'::
        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})
    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).
    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_support : bool
            Whether to display confidence values, if present on the tree.
        ax : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_support`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
    """

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MEDICCPlotError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    if normal_name is not None and len(list(input_tree.find_clades(normal_name))) == 0:
        raise MEDICCPlotError(f'Normal sample "{normal_name}" was not found in tree')

    import matplotlib.collections as mpcollections
    if ax is None:
        nsamp = len(list(input_tree.find_clades()))
        plot_height = height_scale * nsamp * 0.25
        max_leaf_to_root_distances = np.max([np.sum([x.branch_length for x in input_tree.get_path(leaf)])
                            for leaf in input_tree.get_terminals()])
        plot_width = 5 + np.max([0, width_scale * np.log10(max_leaf_to_root_distances / 100) * 5])

        # maximum figure size is 250x250 inches
        fig, ax = plt.subplots(figsize=(min(250, plot_width), min(250, plot_height)))

    label_func=label_func if label_func is not None else lambda x: x

    # options for displaying label colors.
    if label_colors is not None:
        if callable(label_colors):
            def get_label_color(label):
                return label_colors(label)
        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")
    else:
        clade_colors = {}
        for sample in [x.name for x in list(input_tree.find_clades(''))]:
            ## determine if sample is terminal
            is_terminal = True
            matches = list(input_tree.find_clades(sample))
            if len(matches) > 0:
                clade = matches[0]
                is_terminal = clade.is_terminal()
            ## determine if sample is normal
            clade_colors[sample] = COL_MARKER_TERMINAL
            if not is_terminal:
                clade_colors[sample] = COL_MARKER_INTERNAL
            if sample == normal_name:
                clade_colors[sample] = COL_MARKER_NORMAL
        
        def get_label_color(label):
            return clade_colors.get(label, "black")

    if marker_size is None:
        marker_size = TREE_MARKER_SIZE
    marker_func=lambda x: (marker_size, get_label_color(x.name)) if x.name is not None else None

    ax.axes.get_yaxis().set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True, prune=None))
    ax.xaxis.set_tick_params(labelsize=XLABEL_TICK_SIZE)
    ax.xaxis.label.set_size(XLABEL_FONT_SIZE)
    ax.set_title(title, x=0.01, y=1.0, ha='left', va='bottom',
                fontweight='bold', fontsize=16, zorder=10)
    x_posns = _get_x_positions(input_tree)
    y_posns = _get_y_positions(input_tree, adjust=not hide_internal_nodes, normal_name=normal_name)

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def value_to_str(value):
        if value is None or value == 0:
            return None
        elif int(value) == value:
            return str(int(value))
        else:
            return str(value)

    if not branch_labels:
        if show_branch_lengths:
            def format_branch_label(x): 
                return value_to_str(np.round(x.branch_length, 1)) if x.name != 'root' and x.name is not None else None
        else:
            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        def format_branch_label(clade):
            return value_to_str(branch_labels(clade))

    if show_branch_support:
        def format_support_value(clade):
            if clade.name == 'root' or clade.name is None:
                return None
            try:
                confidences = clade.confidences
            # phyloXML supports multiple confidences
            except AttributeError:
                pass
            else:
                return "/".join(value_to_str(cnf.value) for cnf in confidences)
            if clade.confidence is not None:
                return value_to_str(clade.confidence)
            return None


    def draw_clade_lines(
        use_linecollection=False,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            ax.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif not use_linecollection and orientation == "vertical":
            ax.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add node marker
        if marker_func is not None:
            marker = marker_func(clade)
            if marker is not None and clade is not None and not(hide_internal_nodes and not clade.is_terminal()):
                marker_size, marker_col = marker_func(clade)
                ax.scatter(x_here, y_here, s=marker_size, c=marker_col, zorder=3)

        # Add node/taxon labels
        label = label_func(str(clade.name))
        ax_scale = ax.get_xlim()[1] - ax.get_xlim()[0]

        if label not in (None, clade.__class__.__name__) and \
                not (hide_internal_nodes and not clade.is_terminal()):
            ax.text(
                x_here + min(0.02*ax_scale, 1),
                y_here,
                f" {label}",
                verticalalignment="center",
                color=get_label_color(label),
            )
        # Add label above the branch
        conf_label = format_branch_label(clade)
        if conf_label:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                conf_label,
                fontsize="small",
                horizontalalignment="center",
            )
        # Add support below the branch
        if show_branch_support:
            support_value = format_support_value(clade)
            if support_value:
                ax.text(
                    0.5 * (x_start + x_here),
                    y_here + 0.25,
                    support_value + '%',
                    fontsize="small",
                    color='grey',
                    horizontalalignment="center",
                )
        # Add Events list
        if show_events and clade.events is not None:
            ax.text(
                0.5 * (x_start + x_here),
                y_here - 0.15,
                clade.events,
                fontsize="small",
                color=COL_MARKER_NORMAL,
                horizontalalignment="center",
            )
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                use_linecollection=True,
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    if line_width is None:
        line_width = plt.rcParams["lines.linewidth"]
    draw_clade(input_tree.root, 0, "k", line_width)

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        ax.add_collection(i)
    for i in vertical_linecollections:
        ax.add_collection(i)

    ax.set_xlabel("branch length")
    ax.set_ylabel("taxa")

    # Add margins around the `tree` to prevent overlapping the ax
    xmax = max(x_posns.values())
    #ax.set_xlim(-0.05 * xmax, 1.25 * xmax)
    ax.set_xlim(-0.05 * xmax, 1.05 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ax.set_ylim(max(y_posns.values()) + 0.5, 0.5)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                f'Keyword argument "{key}={value}" is not in the format '
                'pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),'
                ' or pyplot_option_name=(dict) '
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if output_name is not None:
        plt.savefig(output_name + ".png", bbox_inches='tight')

    return plt.gcf()


def plot_cn_heatmap(input_df, final_tree=None, y_posns=None, cmax=None, total_copy_numbers=False,
                    alleles=['cn_a', 'cn_b'], tree_width_ratio=1, cbar_width_ratio=0.05, figsize=(20, 10),
                    tree_line_width=0.5, tree_marker_size=0, show_internal_nodes=False, show_branch_support=False, title='',
                    tree_label_colors=None, tree_label_func=None, cmap='coolwarm', normal_name='diploid',
                    ignore_segment_lengths=False):

    input_df = input_df[alleles].copy()
    # if len(np.intersect1d(cur_sample_labels, input_df.index.get_level_values('sample_id').unique())) != len(cur_sample_labels):
    #     raise MEDICCPlotError("tree nodes and labels in dataframe are not the same")

    if not isinstance(alleles, list) and not isinstance(alleles, tuple):
        alleles = [alleles]
    nr_alleles = len(alleles)

    if cmax is None:
        cmax = np.max(input_df[alleles].values.astype(int))

    if final_tree is None:
        fig, axs = plt.subplots(figsize=figsize, ncols=1+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': nr_alleles*[1] + [cbar_width_ratio]})

        cur_sample_labels = (input_df.index.get_level_values('sample_id').unique())
        if not show_internal_nodes:
            logger.warning('No tree provided, so "show_internal_nodes=False" is ignored')
        if y_posns is None:
            y_posns = {s: i for i, s in enumerate(cur_sample_labels)}

        cn_axes = axs[:-1]
        cn_axes[0].set_title(title, x=0, y=1, ha='left', va='bottom', pad=20,
                             fontweight='bold', fontsize=16, zorder=10)
    else:
        fig, axs = plt.subplots(figsize=figsize, ncols=2+nr_alleles, sharey=False,
                                gridspec_kw={'width_ratios': [tree_width_ratio] + nr_alleles*[1] + [cbar_width_ratio]})
        tree_ax = axs[0]
        cn_axes = axs[1:-1]

        if show_internal_nodes:
            cur_sample_labels = np.array([x.name for x in list(final_tree.find_clades()) if x.name is not None])
        else:
            cur_sample_labels = np.array([x.name for x in final_tree.get_terminals()])

        y_posns = {k.name:v for k, v in _get_y_positions(final_tree, adjust=show_internal_nodes, normal_name=normal_name).items()}
        
        _ = plot_tree(final_tree, ax=tree_ax, normal_name=normal_name,
                      label_func=tree_label_func if tree_label_func is not None else lambda x: '',
                      hide_internal_nodes=(not show_internal_nodes), show_branch_lengths=False, 
                      show_events=False, show_branch_support=show_branch_support,
                      line_width=tree_line_width, marker_size=tree_marker_size,
                      title=title, label_colors=tree_label_colors)
        tree_ax.set_axis_off()
        tree_ax.set_axis_off()
        fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0.0, wspace=100)
    cax = axs[-1]

    gaps = (input_df.loc[cur_sample_labels[0]].eval('start') - 
            np.roll(input_df.loc[cur_sample_labels[0]].eval('end'), 1)).values
    total_gaps = gaps[gaps>0].sum()
    if total_gaps > 1e8:
        logger.warning(f"Total of {total_gaps:.1e} bp gaps in the segmentation. These missing "
                    "segments are not reflected in the plot!")

    ind = [y_posns.get(x, -1) for x in cur_sample_labels]

    cur_sample_labels = cur_sample_labels[np.argsort(ind)]
    vcenter = 2 if total_copy_numbers else 1
    color_norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=vcenter, vmax=cmax if cmax > vcenter else vcenter+1)

    chr_ends = input_df.loc[cur_sample_labels[0]].copy()
    chr_ends['end_pos'] = np.cumsum([1]*len(chr_ends))
    chr_ends = chr_ends.reset_index().groupby('chrom').max()['end_pos']
    chr_ends = chr_ends.dropna()
    chr_ends = chr_ends.astype(int)


    if ignore_segment_lengths:
        x_pos = np.arange(len(input_df.loc[cur_sample_labels].astype(int).unstack('sample_id'))+1)
    else:
        x_pos = np.append([0], np.cumsum(input_df.loc[cur_sample_labels].astype(int).unstack(
            'sample_id').loc[:, (alleles[0])].loc[:, cur_sample_labels].eval('end+1-start').values))
    y_pos = np.arange(len(cur_sample_labels)+1)+0.5

    for ax, allele in zip(cn_axes, alleles):
        im = ax.pcolormesh(x_pos, y_pos,
                        input_df.loc[cur_sample_labels].astype(int).unstack(
                            'sample_id').loc[:, (allele)].loc[:, cur_sample_labels].values.T,
                        cmap=cmap,
                        norm=color_norm)

        for _, line in chr_ends.items():
            ax.axvline(x_pos[line], color='black', linewidth=0.75)
        
        xtick_pos = np.append([0], x_pos[chr_ends.values][:-1])
        xtick_pos = (xtick_pos + np.roll(xtick_pos, -1))/2
        xtick_pos[-1] += x_pos[-1]/2
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels([x[3:] for x in chr_ends.index], ha='center', rotation=90, va='bottom')
        ax.tick_params(width=0)
        ax.xaxis.set_tick_params(labelbottom=False, labeltop=True, bottom=False)
        ax.set_yticks([])

    cax.pcolormesh([0, 1],
                   np.arange(0, cmax+2),
                   np.arange(0, cmax+1)[:, np.newaxis],
                   cmap=cmap,
                   norm=color_norm)

    cax.set_xticks([])
    cax.set_yticks(np.arange(0, cmax+1)+0.5)
    if np.max(input_df.values.astype(int)) > cmax:
        cax.set_yticklabels([str(x) + '+' if x == cmax else str(x)
                             for x in np.arange(0, cmax+1)], ha='left')
    else:
        cax.set_yticklabels(np.arange(0, cmax+1), ha='left')
    cax.yaxis.set_tick_params(left=False, labelleft=False, labelright=True)

    for ax in axs[:-1]:
        ax.set_ylim(len(cur_sample_labels)+0.5, 0.5)

    return fig


def compute_cn_change(df, tree, normal_name='diploid'):
    """Compute the copy-number changes per segment in all branches

    Args:
        df (pandas.DataFrame): DataFrame containing the copy-numbers of samples and internal nodes
        tree (Bio.Phylo.Tree): Phylogenetic tree
        normal_name (str, optional): Name of normal sample. Defaults to 'diploid'.

    Returns:
        pandas.DataFrame: DataFrame containing the copy-number changes
    """    
    cn_change = df.copy()
    alleles = cn_change.columns
    for allele in alleles:
        cn_change[allele] = cn_change[allele].astype('int')

    clades = [clade for clade in tree.find_clades(order = "postorder") if clade.name is not None and clade.name != normal_name]
    for clade in clades:
        for child in clade.clades:
            cn_change.loc[child.name, alleles] = cn_change.loc[child.name, alleles].values - cn_change.loc[clade.name, alleles].values
    cn_change.loc[clades[-1].name, alleles] = cn_change.loc[clades[-1].name, alleles].values - cn_change.loc[normal_name, alleles].values
    cn_change.loc[normal_name, alleles] = 0

    return cn_change

def plot_nni_trace(trace):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(trace, '-o', alpha=0.7, markersize=4)
    ax.set_xlabel('NNI Sweep')
    ax.set_ylabel('Sum of Branch Lengths')
    ax.set_title('NNI Trace')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    return fig

class MEDICCPlotError(Exception):
    pass
