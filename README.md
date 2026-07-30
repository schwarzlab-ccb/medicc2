# MEDICC2 - Whole-genome doubling-aware copy number phylogenies for cancer evolution

[![PyPI](https://img.shields.io/pypi/v/medicc2?color=green)](https://pypi.org/project/medicc2/)
[![Conda](https://img.shields.io/conda/v/bioconda/medicc2?color=green)](https://anaconda.org/bioconda/medicc2)
[![Tests](https://github.com/ICCB-Cologne/medicc2/actions/workflows/tests.yml/badge.svg)](https://github.com/ICCB-Cologne/medicc2/actions/workflows/tests.yml)

For more information see the accompanying publication [Whole-genome doubling-aware copy number phylogenies for cancer evolution with MEDICC2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02794-9).

# Installation
Install MEDICC2 via conda (recommended), pip or from source. MEDICC2 was developed and tested on unix-built systems (Linux and MacOS). For Windows users we recommended WSL2.

Note that the notebooks and examples are not included when installing from conda or pip.

For all installation methods you need to make sure to have a working version of the GNU Compiler Collection (`gcc`, `gxx` as well as related packages such as `libgcc-ng`) installed. Note that MEDICC2 requires Cython version 3.0 or later.


## Installation via conda (recommended)
It is best to use a dedicated conda environment for your MEDICC2 installation with `conda create -n medicc_env`.

After activating the environment with `conda activate medicc_env` you can install MEDICC2 via `conda install -c bioconda -c conda-forge medicc2`.

## Installation via pip
As MEDICC2 relies on OpenFST version 1.8.2 which is not packaged on PyPi you have to first install it using conda with `conda install -c conda-forge openfst=1.8.4`. Next you can install MEDICC2 via `pip install medicc2`.

## Installation from source
Clone the MEDICC2 repository and its submodules using `git clone --recursive https://github.com/ICCB-Cologne/medicc2.git`. It is important to use the `--recursive` flag to also download the modified OpenFST submodule.

All dependencies including OpenFST (v1.8.2) should be directly installable via conda. A yaml file with a suggested MEDICC2 conda environment is provided in 'doc/medicc2.yml'. You can create a new conda environment with all requirements using `conda env create -f doc/medicc2.yml -n medicc_env`.

Then, inside the `medicc2` folder, run `pip install .` to install MEDICC2 to your environment.

## Development mode
If you want to make changes to the MEDICC2 source code, install MEDICC2 as explained above in "Installation from source" but install in editable mode (with `pip install -e .`). Any changes to the local files will now be reflected in your environment.

# Usage
After installing MEDICC2, you can use MEDICC2 functions in python scripts (through `import medicc`) and from the command line. General usage from the command line is `medicc2 path/to/input/file path/to/output/folder`. Run `medicc2 --help` for information on optional arguments.

Logging settings can be changed using the `medicc/logging_conf.yaml` file with the standard python logging syntax.

## Command line Flags

* `input_file`: path to the input file
* `output_dir`: path to the output folder
* `--version`: Print version information
* `--input-type`, `-i`: Choose the type of input: f for FASTA, t for TSV. Default: 'TSV'
* `--input-allele-columns`, `-a`: Name of the CN columns (comma separated) if using TSV input format. This also adjusts the number of alleles considered (min. 1, max. 2). Default: 'cn_a, cn_b'
* `--input-chr-separator`: Character used to separate chromosomes in the input data (condensed FASTA only). Default: 'X'
* `--tree`: Do not reconstruct tree, use provided tree instead (in newick format) and only perform ancestral reconstruction. Default: None
* `--topology-only`, `-s`: Output only tree topology, without reconstructing ancestors. Default: False
* `--length-aware-ancestors`, `-la`: Use the event-counting FST for the upper pass and the length-encoding FST for lower-pass ancestral selection, while keeping the MEDICC distance matrix and reported branch lengths event-counted. Default: False
* `--nni-mode`: Refine the tree topology with Nearest Neighbor Interchange (NNI) search after the initial neighbor-joining/supplied tree. See "Nearest Neighbor Interchange (NNI) search" below. Default: False
* `--nni-max-topology`: Hard cap on the total number of NNI neighbor evaluations across all search sweeps (an evaluation budget, not a count of distinct topologies). Only relevant with `--nni-mode`. Default: 20000
* `--nni-export-all-topology`: Export every co-optimal topology NNI found, instead of just the first one, each into its own output subfolder. Only relevant with `--nni-mode`; cannot be combined with bootstrapping. Default: False
* `--nni-complete-log`: Write a complete per-sweep trace of NNI mode (topology and score at every step) to `<output_dir>/nni-trace/`. Only relevant with `--nni-mode`. Default: False
* `--normal-name`, `-n`: ID of the sample to be treated as the normal sample. Trees are rooted at this sample for ancestral reconstruction. If the sample ID is not found, an artificial normal sample of the same name is created with CN states = 1 for each allele. Default: 'diploid'
* `--exclude-samples`, `-x`: Comma separated list of sample IDs to exclude. Default: None
* `--filter-segment-length`: Removes segments that are smaller than specified length (measured in bp's). Default: None
* `--bootstrap-method`: Bootstrap method. Has to be either 'chr-wise' or 'segment-wise'. Default: 'chr-wise'
* `--bootstrap-nr`: Number of bootstrap runs to perform. Default: None
* `--prefix`, '-p': Output prefix to be used. None uses input filename. Default: None
* `--no-wgd`: Disable whole-genome doubling events. Default: False
* `--plot`: Type of copy-number plot to save. 'bars' is recommended for <50 samples, heatmap for more samples, 'auto' will decide based on the number of samples, 'both' will plot both and 'none' will plot neither. (default: auto).
* `--no-plot-tree`: Disable plotting of tree figures. Default: False
* `--total-copy-numbers`: Run for total copy number data instead of allele-specific data. Default: False
* `-j`, `--n-cores`: Number of cores to run on. Default: None
* `--events`: Whether to infer copy-number events. See section "Event Reconstruction" below
* `--chromosomes-bed`: BED file for chromosome regions to compare copy-number events to
* `--regions-bed`: BED file for regions of interest to compare copy-number events to
* `-v`, `--verbose`: Enable verbose output. Default: False
* `-vv`, `--debug`: Enable more verbose output Default: False
* `--silent`: Hide all output. Default: False
* `--maxcn`: Expert option: maximum CN at which the input is capped. Does not change FST. The maximum possible value is 8. Default: 8
* `--prune-weight`: Expert option: Prune weight in ancestor reconstruction. Values >0 might result in more accurate ancestors but will require more time and memory. Default: 0
* `--fst`: Expert option: path to an alternative FST. **Deprecated** — use `--upper-fst`/`--lower-fst` instead. Default: None
* `--upper-fst`: Expert option: path to an alternative upper FST. Must be given together with `--lower-fst`. Default: None
* `--lower-fst`: Expert option: path to an alternative lower FST. Must be given together with `--upper-fst`. Default: None
* `--fst-chr-separator`: Expert option: character used to separate chromosomes in the FST. Default: 'X'
* `--wgd-x2`: Expert option: Treat WGD as a x2 operation. Default: False
* `--lower-fst-nni`, `-l-nni`: Expert option: also use the lower FST (rather than only the upper FST) while scoring candidates during NNI search itself, not just for the final ancestral reconstruction. Only relevant with `--nni-mode`. Default: False


## Input files
Input files can be either in fasta or tsv format:

* **tsv (recommended):** Files should have the following columns: `sample_id`, `chrom`, `start`, `end` as well as columns for the copy numbers. MEDICC expects the copy number columns to be called `cn_a` and `cn_b`. Using the flag `--input-allele-columns` you can set your own copy number columns. If you want to use total copy numbers, make sure to use the flag `--total-copy-numbers`. Important: MEDICC2 does not create total copy numbers for you. You will have to calculate total copy numbers yourself and then specify the column using the `--input-allele-columns` flag.
* **fasta:** A description file should be provided to MEDICC. This file should include one line per file with the name of the chromosome and the corresponding file names. If fasta files are provided you have to use the flag `--input-type fasta`.

Note that MEDICC2 needs at least 2 non-diploid samples to create a tree.

MEDICC2 follows the BED convention for segment coordinates, i.e. segment start is at 0 and the segment end is non-inclusive.

The folder `examples/simple_example` contains a simple example input both in fasta and tsv format.
The folder `examples/OV03-04` contains a larger example consisting of multiple fasta files. If you want to run MEDICC on this data run `medicc2 examples/OV03-04/OV03-04_descr.txt path/to/output/folder --input-type fasta`.

**Note that MEDICC2 requires a consistent segmentation across all samples**. That means that all segments must have the exact same segments.


## Output files
MEDICC creates the following output files:

* `_final_tree.new`, `_final_tree.xml`, `_final_tree.png`: The final phylogenetic tree in Newick and XML format as well as an image
* `_pairwise_distances.tsv`: A NxN matrix (N being the number of samples) of pairwise distances calculated with the symmetric MEDICC2 distance (from every pair of samples). This matrix is used to create the phylogenetic tree through neighbor joining. Note that this distance matrix is different from the distances in the final tree due to the inferred ancestral states present in tree; to get the pairwise distances from the final tree, we recommend using an external tool such as the "cophenetic.phylo" function from the "ape" package in R.
* `_final_cn_profiles.tsv`: Copy-number profiles of the input as well as the newly internal nodes. Also includes additional information such as whether a gain or loss has happened.
* `_cn_profiles.pdf`: Combined plot of the phylogenetic tree as well as the copy-number profiles of all samples (including the internal nodes)
* `_branch_lengths.tsv`: List of all branches and their corresponding lenghts of the final tree
* `_summary.tsv`: Contains summary information about the created tree. If the `--events` flag was set, this includes the WGD status.

*optional (see "Event Reconstruction" below)*

* `_copynumber_events_df.tsv`: List of all copy-number events detected. The entries for WGD events have non-meaningful values for chrom, cn_child, etc. Note that the events derived are not unambiguous.
* `_events_overlap.tsv`: Overlap of copy-number events with regions of interest

*optional (see "Nearest Neighbor Interchange (NNI) search" below)*

* `nni_trace.png`: Written whenever `--nni-mode` is used. Plots the best score found at each NNI search sweep.
* `nni-trace/nni_trace.tsv`: Written only if `--nni-complete-log` is also set. A complete per-step trace of every topology and score considered during the search.

### Output files with `--nni-export-all-topology`
When `--nni-export-all-topology` is combined with `--nni-mode`, the output directory structure changes: instead of a single set of `_final_tree.*`/`_branch_lengths.tsv`/`_final_cn_profiles.tsv`/`_summary.tsv`/plot files directly in the output directory, MEDICC2 creates one subfolder per co-optimal topology found, numbered `topology_1/`, `topology_2/`, and so on. Each subfolder contains that topology's own tree, branch lengths, copy-number profiles, plots, summary, and (if `--events` was set) copy-number events and regions-of-interest overlap files, using the same filenames as above with `_topology_N` added (e.g. `_final_tree_topology_2.new`, `_final_cn_profiles_topology_2.tsv`) — except `_summary.tsv`, `_copynumber_events_df.tsv` and `_events_overlap.tsv`, which are written per-topology inside each subfolder without the `_topology_N` suffix. Only `_pairwise_distances.tsv` and `nni_trace.png` are shared across all topologies and stay at the top level of the output directory.


## Output plots
Apart from the file `_tree.pdf` which contains the inferred phylogeny, the main plot created by MEDICC is the copy-number plots named either `_cn_profiles.pdf` or `_cn_profiles_heatmap.pdf`.
The left part consists of the inferred phylogenetic tree including the number of events in the branches. The right part is made up of the copy-number profiles of the samples (and potentially the reconstructed ancestral nodes).

There are two kinds of copy-number plots: the bars and the heatmap version. The bars version is most suitable for fewer samples (<50) as more details are visible while the heatmap version is most suitable many samples expected for example in single-cell experiments.
You can toggle the kind of plot MEDICC2 creates with the `--plot` flag (see above).

### Example bars copy-number plot
Example from patient PTX011 from the Gundem et al. Nature 2015. The data can be found in `example/gundem_et_al_2015/`.

![copy-number bars plot for PTX011 Gundem 2015](doc/MEDICC2_cn_plot_example_bars.png)

**Legend**

![legend of copy-number plot](doc/MEDICC2_cn_plot_legend.png)

### Example heatmap copy-number plot
Example from patient PTX011 from the Gundem et al. Nature 2015. The data can be found in `example/gundem_et_al_2015/`.

![copy-number heatmap plot for PTX011 Gundem 2015](doc/MEDICC2_cn_plot_example_heatmap.png)


## Usage examples
For first time users we recommend to have a look at `examples/simple_example` to get an idea of how input data should look like. Then run `medicc2 examples/simple_example/simple_example.tsv path/to/output/folder` as an example of a standard MEDICC run. Finally, the notebook `notebooks/example_workflows.py` shows how the individual functions in the workflow are used.

The notebook `notebooks/bootstrap_demo.py` demonstrates how to use the bootstrapping routine and `notebooks/plot_demo.py` shows how to use the main plotting functions.


## Event Reconstruction
MEDICC2 can create a list of copy-number events in the file `_copynumber_events_df.tsv` which are also displayed in the final copy-number barplot.
These are disabled by default and are enabled using the `--events` flag.

Note, that the inferred events are not unambigous but just one possible solution. In some cases there are multiple possible paths that can result in the same final copy-number state in the same number of steps. Without additional information, MEDICC2 cannot determine which possible path is the right one and thus opts for a path that creates the longest consecutive gains. 
Even though the events inferred by MEDICC2 are not unambigous they are minimal (as in there are no solutions with fewer number of steps) and deterministic (as in multiple runs of MEDICC2 will always return the same events).

Minimal example: *111 -> 232* which can be explained by *gain-gain-gain* + *x-gain-x* or *gain-gain-x* + *x-gain-gain*. MEDICC2 would select the first option.

Therefore subsequent analysis on the reconstructed events can potentially be influenced by artifacts and most be treated with caution.


### Regions of interest
MEDICCC2 can overlap the inferred events with regions of interest such as chromosome arms or oncogenes. 
This process requires the installation of `pyranges` which might be incompatible with newer version of python and/or numpy.
The overlap is turned off by default. You can turn on the overlapping with the ``--chromosomes-bed` and `--regions-bed` flag by providing bed-files with regions of interest. By default MEDICC2 uses hg38 chromosome-arms and a list of genes taken from Davoli et al. Cell 2013. This data is present as BED files in the `medicc/objects` folder. Invoke these using the flags `--chromosomes-bed "default"` and/or `--regions-bed "default"`.
Users can specify regions of interest of their own in BED format by providing the `--chromosomes-bed` or `--regions-bed` flags.


## Nearest Neighbor Interchange (NNI) search
By default, MEDICC2 builds its tree topology with neighbor joining and does not search further. Passing `--nni-mode` additionally refines that topology (or a topology supplied via `--tree`) with a Nearest Neighbor Interchange (NNI) local search: each search step tries swapping neighboring branches and keeps whichever resulting topology has the best (lowest) score, repeating until no swap improves the score further.

* `--nni-max-topology` sets a hard cap on the total number of NNI neighbor evaluations across the whole search (a compute budget, not a count of distinct topologies), in case the search has not converged yet.
* The search can end in a tie between multiple equally good topologies. By default MEDICC2 reports only one of them; `--nni-export-all-topology` exports every co-optimal topology instead, each written to its own `topology_N/` subfolder under the output directory. This cannot be combined with bootstrapping (`--bootstrap-nr`/`--bootstrap-method`), since bootstrapping requires a single reference tree.
* `--nni-complete-log` writes a full per-sweep trace (topology and score at every step) to `<output_dir>/nni-trace/`.
* `--lower-fst-nni`/`-l-nni` and `--length-aware-ancestors`/`-la` are expert options controlling which FST is used during NNI search and ancestral reconstruction respectively; see the flag descriptions above.

NNI mode always starts from whichever tree MEDICC2 has already built (via neighbor joining, or a tree supplied with `--tree`) — there is no separate flag to give NNI a different starting tree.


## Single sample WGD detection
If you are interested in the WGD status of individual samples in your data, have a look at the notebook `notebooks/WGD_detection.ipynb`. The notebooks also includes code to detect the presence of multple WGDs in your samples. By replacing the input data with your data you can easily calculate the WGD status of any copy-number input.


# Issues
If you experience problems with MEDICC2 please [file an issue directly on Bitbucket](https://bitbucket.org/schwarzlab/medicc2/issues/new) or [contact us directly](chenxi.nie@iccb-cologne.org). 

## Known Issues

**Phasing**
We do recommend to phase your input data before using MEDICC2. MEDICC2's own phasing algorithm is only to be used when looking at single samples and should not be used in the case of multiple samples! Here we recommend using Refphase ([Bitbucket](https://bitbucket.org/schwarzlab/refphase/)).
If phasing is not possible for you, working on major/minor configuration works reasonably well in most cases.
For very noisy data, where accurate phasing cannot be guaranteed, you can also try to create total copy numbers and run MEDICC2 with the `--total-copy-numbers` flag.

**long runtime**
MEDICC2 tries to solve an NP-hard problem by inferring a symmetric distance between samples and therefore has a higher runtime than other tools than only compute the asymmetric (and less accurate) distance between samples.

In order to improve runtime, you should first run MEDICC2 with multiple cores using the `--n-cores` flag. Using multiple cores will roughly lead to a decrease in runtime linear w.r.t. numbers of cores used (depending on your system architecture).

Next, you can remove duplicate and diploid cells. Especially in the case of 100s to 1000s single-cell samples, there are oftentimes multiple copies of cells with the same copy-number profiles as well as cells that are (almost) purely diploid. Removing those will (in most cases) not alter the results but decrease runtime.

Finally, MEDICC2's runtime scales with the number of segments used. If you use copy-number bins, try to increase the bin-size for decreased runtime as well as merge neighboring bins that are equal in copy-number states across all samples. (For example a chromosome that is purely diploid in all samples should be represented as a single segment rather than multiple bins).
If you instead create a minimum consistent segmentation, be aware of individual samples with many breakpoints that will drive up your total number of breakpoints and therefore number of segments. If individual samples have an excessive amount of breakpoints, it is best to remove them before creating a minimum consistent segmentation.

**Noisy segments**
Small faulty or noisy segments can have a strong effect on the distances MEDICC2 calculates between samples and therefore the resulting tree.
This is because MEDICC2 counts all segments equally in order appropriatlely take focal events into account. 
If the resulting and the inferred events look strange to you, you can replot the tree and copy-number profiles using the function `plot_cn_profiles` setting `ignore_segment_lengths=True` (see the notebook `notebooks/plot_demo.py` for usage examples) in order to investigate small segments that might not have been visible in the original plot.
If you are unsure about the copy-number profiles we recommened to filter small segments.

**Taxon imbalance**
If your data contains 100s to 1000s samples with a few distinct subgroups, an imbalance in the number of samples per subgroups might lead to an incorrect tree (e.g. 50 samples of subclone A and 1000 samples each of subclone B and C).
This is a known problem in phylogeny called *taxon imbalance* or *taxon sampling*. If you have multiple, clearly separable subgroups in your data we recommoned either subsampling over-represented groups or upsampling under-represented groups to gauge the effect of taxon imbalance.

**Running out of memory / bad_alloc error**
If MEDICC2 terminates with the following error `terminate called after throwing an instance of 'std::bad_alloc'` or your machine runs out of memory this hints towards an issue with the FST.
Rerun MEDICC2 with the `-vv` flag to enable extended logging. If the error occurs during the ancestral reconstruction routine, the issue is related to OpenFST which is the FST library employed by MEDICC2 and cannot be easily solved by us.
This issue can be related to small bin sizes (and therefore a large number of segments). Increasing the binsize (although decreasing accuracy) solves this issue most of the time.
You can also try to remove the sample that led to the error (see the extended logs for this). 

**The output plots are not like I expected**
Maybe you need to set the `--plot` flag. By default, `--plot` is set to auto which means that it plots different figures depending on the number of samples in the data (threshold is 50); see above.

**Faulty event reconstruction**
Sometimes MEDICC2 will pass out the following warning: *Event recreation was faulty*. This means that the events in the
`_cn_events_df.tsv` file will not be accurate. If you selected total copy number this will mainly be due to multiple WGDs
in a single node. Please get in contact with us if the problem prevails even without the `--total-copy-numbers` flag.

**missing segments / gaps in the segmentation**
MEDICC2 will assume that the segmentation is gap-less, i.e. that gaps between neighboring segments are neglible. If your data contains large gaps this might affect the performance of MEDICC2 as it might incorrectly jointly mutate segments that are actually separated.


# Bugs, feature requests and contact
You can report bugs and request features directly in [Github](https://github.com/schwarzlab-ccb/medicc2/issues) or contact us via at *chenxi.nie@iccb-cologne.org*.


# License
MEDICC2 is available under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html). It contains modified code of the *pywrapfst* Python module from [OpenFST](http://www.openfst.org/) as permitted by the [Apache 2](http://www.apache.org/licenses/LICENSE-2.0) license.


# Kaufmann et al. 2022 Publication
The MEDICC2 model has been published in 2022 with Genome Biology: [MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution](https://doi.org/10.1186/s13059-022-02794-9).

## Figures
Both the required data and the scripts to create all Figures from the publication are stored in the commit version *9b400ef* of the MEDICC2 repository available on [Zenodo](https://zenodo.org/record/7300106) in the folder `Figures_Kaufmann_et_al_2021/`.

## Simulation Validation
The simulated copy-number profiles used in the MEDICC2 publication were created using our simulation framework [Simphyni](https://bitbucket.org/schwarzlab/simphyni/src/main/). The Simphyni repository contains a notebook to recreate the exact data used in the publication.

## Please cite
Kaufmann, T.L., Petkovic, M., Watkins, T.B.K. et al.  
**MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution**.  
Genome Biol 23, 241 (2022). https://doi.org/10.1186/s13059-022-02794-9

Schwarz RF, Trinh A, Sipos B, Brenton JD, Goldman N, Markowetz F.  
**Phylogenetic quantification of intra-tumour heterogeneity.**  
PLoS Comput Biol. 2014 Apr 17;10(4):e1003535. doi: 10.1371/journal.pcbi.1003535.
