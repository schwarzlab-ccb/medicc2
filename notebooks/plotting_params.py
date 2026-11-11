import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation 
import seaborn as sns

plotting_params = {
    'WIDTH_FULL': 12,
    'WIDTH_HALF': 6,
    'ASPECT_RATIO': 4/3,
    'FONTSIZE_LARGE': 12,
    'FONTSIZE_MEDIUM': 10,
    'FONTSIZE_SMALL': 10, # was 8
    'FONTSIZE_TINY': 8, # was 5
    'LINEWIDTH': 3,
    'MARKERSIZE_SMALL': 3,
    'MARKERSIZE_MEDIUM': 5,
    'MARKERSIZE_LARGE': 10,
    'LINEWIDTH_SMALL': 1}

color_palette = sns.color_palette([
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf"
])


def set_plotting_params():
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=plotting_params['FONTSIZE_MEDIUM'])
    plt.rc('axes', titlesize=plotting_params['FONTSIZE_LARGE'])
    plt.rc('axes', labelsize=plotting_params['FONTSIZE_MEDIUM'])
    plt.rc('xtick', labelsize=plotting_params['FONTSIZE_SMALL'])
    plt.rc('ytick', labelsize=plotting_params['FONTSIZE_SMALL'])
    plt.rc('legend', fontsize=plotting_params['FONTSIZE_SMALL'])
    plt.rc('legend', frameon=False)
    plt.rc('figure', titlesize=plotting_params['FONTSIZE_LARGE'])

    sns.set_palette(sns.color_palette(color_palette))

    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False


def label_axes(axs):
    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
              
    for ax, l in zip(axs.flat, letters):
        trans = ScaledTranslation(-20/72, 7/72, plt.gcf().dpi_scale_trans)
        ax.text(-0.025, 1.0, l, transform=ax.transAxes + trans,
            # fontsize=30, color='black', va='bottom', ha='left', fontfamily='sans-serif', fontweight='bold')
            fontsize=30, color='black', va='top', ha='right', fontfamily='sans-serif', fontweight='bold')
