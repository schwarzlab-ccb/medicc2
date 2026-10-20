"""MEDICC reference CLI for the full L_LOH/two-WGD model."""

import argparse
import csv
import gzip
import json
import logging
import math
from pathlib import Path
import time
import medicc
from Bio import Phylo
from .core import reconstruct


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--version",
        action="version",
        version=medicc.__version__ + " full-fst-l-loh-wgd-2",
    )
    parser.add_argument("--normal-name", "-n", default="diploid")
    parser.add_argument(
        "--input-type", "-i", choices=["t", "tsv", "f", "fasta"], default="tsv"
    )
    parser.add_argument("--input-allele-columns", "-a", default="cn_a,cn_b")
    parser.add_argument(
        "--tree", help="Use a supplied Newick topology; skip distance inference"
    )
    parser.add_argument(
        "--nni-mode",
        action="store_true",
        help="Run MEDICC NNI after FastME or from --tree",
    )
    parser.add_argument(
        "--nni-max-topology",
        type=int,
        default=None,
        help="Optional NNI evaluation limit; default unlimited",
    )
    parser.add_argument("--nni-max-frontier", type=int, default=1)
    parser.add_argument(
        "--nni-traverse-plateaus",
        action="store_true",
        help="Default stops after a sweep without strict improvement",
    )
    parser.add_argument("--events", action="store_true",
                        help="Call minimum-cost events using the capped full L_LOH/two-WGD model")
    parser.add_argument("--chromosomes-bed", default="none",
                        help="With --events: chromosome-arm BED file, default (hg38), or none")
    parser.add_argument("--regions-bed", default="none",
                        help="With --events: named-region BED file, default (Davoli genes), or none")
    parser.add_argument("--topology-only", "-s", action="store_true")
    parser.add_argument("--n-cores", "-j", type=int, default=1)
    parser.add_argument("--fastme", default="fastme")
    parser.add_argument("--prefix", "-p")
    parser.add_argument("--maxcn", type=int, default=8, choices=range(1, 9))
    parser.add_argument(
        "--plot", choices=["none", "auto", "bars", "heatmap", "both"], default="none"
    )
    parser.add_argument("--no-plot-tree", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)
    if args.n_cores < 1:
        parser.error("--n-cores must be positive")
    if args.events and args.topology_only:
        parser.error("--events requires ancestors")
    if args.nni_mode and args.topology_only:
        parser.error("--nni-mode requires ancestors")
    if args.nni_max_topology is not None and args.nni_max_topology < 1:
        parser.error("--nni-max-topology must be positive")
    if args.nni_max_frontier < 1:
        parser.error("--nni-max-frontier must be positive")
    if args.topology_only and args.plot != "none":
        parser.error("Use --plot none with --topology-only")
    annotate = args.chromosomes_bed != "none" or args.regions_bed != "none"
    if annotate:
        if not args.events:
            parser.error("BED annotations require --events")
        try:
            import pyranges
        except ImportError:
            parser.error("BED annotations require PyRanges; install medicc2-full-fst[events] "
                         "or run: python -m pip install 'pyranges>=0.1,<0.2'")
        objects = Path(__file__).resolve().parents[1] / "objects"
        for option, default in [("chromosomes_bed", "hg38_chromosome_arms.bed"),
                                ("regions_bed", "Davoli_2013_TSG_OG_genes.bed")]:
            value = getattr(args, option)
            if value == "none":
                continue
            path = objects / default if value == "default" else Path(value).expanduser()
            if not path.is_file():
                parser.error(f"BED file does not exist: {path}")
            setattr(args, option, str(path))
    columns = [s.strip() for s in args.input_allele_columns.split(",")]
    if len(columns) != 2 or len(set(columns)) != 2:
        parser.error("This reference CLI expects two distinct allele columns")
    prefix = args.prefix or Path(args.input_file).stem
    if Path(prefix).name != prefix or prefix in ("", ".", ".."):
        parser.error("Prefix must be a filename component")
    out = args.output_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / (prefix + "_summary.json")).exists():
        parser.error(
            "This output prefix already exists; choose a new output directory or prefix"
        )
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)
    start = time.perf_counter()
    frame = medicc.io.read_and_parse_input_data(
        args.input_file,
        normal_name=args.normal_name,
        input_type=args.input_type,
        allele_columns=columns,
        maxcn=args.maxcn,
    )
    result = reconstruct(
        frame,
        normal_name=args.normal_name,
        input_tree=args.tree,
        n_cores=args.n_cores,
        fastme=args.fastme,
        nni=args.nni_mode,
        topology_only=args.topology_only,
        events=args.events,
        output_dir=out,
        allele_columns=columns,
        nni_max_iter=args.nni_max_topology
        if args.nni_max_topology is not None
        else math.inf,
        nni_max_frontier=args.nni_max_frontier,
        nni_stop_on_plateau=not args.nni_traverse_plateaus,
    )
    if result.pairwise_distances is not None:
        medicc.io.write_pairwise_distances(
            result.sample_labels,
            result.pairwise_distances,
            str(out / (prefix + "_pairwise_distances")),
        )
    Phylo.write(result.initial_tree, out / (prefix + "_initial_tree.new"), "newick")
    medicc.io.write_tree_files(
        result.tree,
        str(out / (prefix + "_final_tree")),
        plot_tree=not args.no_plot_tree,
        normal_name=args.normal_name,
    )
    medicc.io.write_branch_lengths(
        result.tree, str(out / (prefix + "_branch_lengths.tsv"))
    )
    if result.profiles is not None:
        result.copy_numbers.to_csv(out / (prefix + "_final_cn_profiles.tsv"), sep="\t")
        with gzip.open(out / (prefix + "_profiles.json.gz"), "wt") as handle:
            json.dump(result.profiles, handle)
    if result.events is not None:
        result.events.to_csv(out / (prefix + "_copynumber_events_df.tsv"), sep="\t")
    if annotate:
        from medicc.event_reconstruction import overlap_events
        overlaps = overlap_events(events_df=result.events.copy(), alleles=columns,
                                  normal_name=args.normal_name,
                                  chromosome_bed=args.chromosomes_bed,
                                  regions_bed=args.regions_bed)
        overlaps.to_csv(out / (prefix + "_events_overlap.tsv"), sep="\t")
    if result.nni_result:
        with (out / (prefix + "_nni_steps.tsv")).open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["step", "tree", "score"])
            writer.writerows(result.nni_result["step_records"])
    summary = dict(
        result.details,
        elapsed_seconds=time.perf_counter() - start,
        implementation=medicc.__file__,
        status="ok",
    )
    (out / (prefix + "_summary.json")).write_text(json.dumps(summary, indent=2) + "\n")
    if args.plot != "none":
        import matplotlib.pyplot as plt

        styles = ["bars", "heatmap"] if args.plot == "both" else [args.plot]
        if styles == ["auto"]:
            styles = ["bars" if len(result.sample_labels) < 50 else "heatmap"]
        for style in styles:
            if style == "bars":
                figure = medicc.plot.plot_cn_profiles(
                    result.copy_numbers,
                    input_tree=result.tree,
                    normal_name=args.normal_name,
                    allele_columns=columns,
                )
            else:
                figure = medicc.plot.plot_cn_heatmap(
                    result.copy_numbers,
                    final_tree=result.tree,
                    alleles=columns,
                    normal_name=args.normal_name,
                    total_copy_numbers=False,
                )
            figure.savefig(
                out / (prefix + "_cn_" + style + ".pdf"), bbox_inches="tight"
            )
            plt.close(figure)
    print(
        f"MEDICC full-FST finished: {out}; total branch score {summary['tree_length']}"
    )
