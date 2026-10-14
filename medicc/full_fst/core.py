"""MEDICC reconstruction with the full L_LOH model and up to two WGDs per edge.

The native kernel is exact for its fixed model. Ancestral messages retain all
local minima (prune_weight=0), as in the experimental MEDICC reconstruction.
"""

from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import csv
import math
import re
import shutil
import subprocess
import tempfile
import threading
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from .model import Model
from .pair_scorer import Model as DistanceModel
from .scheduler import reconstruct as reconstruct_messages
from . import shortening


@dataclass
class Result:
    sample_labels: list
    pairwise_distances: pd.DataFrame | None
    initial_tree: Tree
    tree: Tree
    profiles: dict | None
    copy_numbers: pd.DataFrame
    details: dict
    nni_result: dict | None = None


def encode(input_df, normal_name, allele_columns=None):
    """Encode the standard MEDICC dataframe without constructing string FSTs."""
    if list(input_df.index.names) != ["sample_id", "chrom", "start", "end"]:
        raise ValueError("Expected MEDICC index: sample_id, chrom, start, end")
    if not input_df.index.is_unique:
        raise ValueError("Duplicate input sample/bin rows")
    columns = list(allele_columns or input_df.columns)
    if len(columns) not in (1, 2):
        raise ValueError("Expected one or two allele columns")
    labels = list(input_df.index.get_level_values("sample_id").unique())
    if any(not isinstance(s, str) or not s for s in labels):
        raise ValueError("Sample IDs must be nonempty strings")
    if normal_name not in labels:
        raise ValueError(
            "Normal sample missing; use medicc.io.read_and_parse_input_data to add it"
        )
    labels = [normal_name] + sorted(
        (s for s in labels if s != normal_name),
        key=lambda s: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)],
    )
    if len(labels) < 3:
        raise ValueError("At least two tumour samples are required")
    coordinates = input_df.loc[normal_name].sort_index().index
    genomes = {}
    for sample in labels:
        frame = input_df.loc[sample].sort_index()
        if not frame.index.equals(coordinates) or len(frame) == 0:
            raise ValueError(f"Inconsistent or empty bin grid: {sample}")
        words = []
        for allele in columns:
            previous = None
            chunk = []
            for (chrom, start, end), value in frame[allele].items():
                word = str(value)
                if word not in "012345678" or len(word) != 1:
                    raise ValueError(
                        f"CN must be an integer from 0 through 8: {sample}: {value}"
                    )
                if end < start:
                    raise ValueError("Segment end precedes start")
                if previous is not None and chrom != previous:
                    words.append("".join(chunk))
                    chunk = []
                chunk.append(word)
                previous = chrom
            words.append("".join(chunk))
        genomes[sample] = "X".join(words)
    return genomes, coordinates, columns


def distance_matrix(genomes, workers=1, progress=None):
    """Symmetric kernel distances; each worker owns its native cache."""
    if workers < 1:
        raise ValueError("workers must be positive")
    progress = progress or (lambda event: None)
    labels = list(genomes)
    n = len(labels)
    values = np.zeros((n, n), dtype=float)
    prepared = {name: shortening.prepare(word) for name, word in genomes.items()}
    local = threading.local()
    models = []
    lock = threading.Lock()

    def score_row(i):
        if not hasattr(local, "model"):
            local.model = DistanceModel()
            with lock:
                models.append(local.model)
        row = []
        columns = 0
        for j in range(i + 1, n):
            x, y = shortening.pair(prepared[labels[i]], prepared[labels[j]])
            value = local.model.score(x, y)
            if not math.isfinite(value):
                raise ValueError(f"Infeasible pair: {labels[i]}, {labels[j]}")
            row.append((j, value))
            columns += len(x)
        local.model.clear()
        return i, row, columns

    pairs = columns = 0
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for future in as_completed(
                [pool.submit(score_row, i) for i in range(n - 1)]
            ):
                i, row, shortened = future.result()
                pairs += len(row)
                columns += shortened
                for j, value in row:
                    values[i, j] = values[j, i] = value
                progress(
                    dict(
                        phase="distances",
                        completed_pairs=pairs,
                        total_pairs=n * (n - 1) // 2,
                    )
                )
    finally:
        for model in models:
            model.close()
    return pd.DataFrame(values, index=labels, columns=labels), dict(
        pairs=pairs,
        mean_shortened_columns=columns / pairs if pairs else 0,
        original_columns=len(next(iter(genomes.values()))) if genomes else 0,
    )


def normalize_tree(tree, normal_name, labels, rename=False):
    """Use the unary normal root expected by MEDICC's incremental NNI search."""
    tree = copy.deepcopy(tree)
    if tree.root.name == normal_name and len(tree.root.clades) == 1:
        root = tree.root.clades[0]
    elif normal_name in {n.name for n in tree.get_terminals()}:
        tree.root_with_outgroup(normal_name)
        rest = [n for n in tree.root.clades if n.name != normal_name]
        if len(rest) != 1:
            raise ValueError("Expected one tumour subtree after rooting on normal")
        root = rest[0]
    else:
        root = tree.root
    if {n.name for n in root.get_terminals()} != set(labels) - {normal_name}:
        raise ValueError("Tree tips do not match observed samples")
    used = {n.name for n in root.find_clades() if n.name is not None} | {normal_name}
    counter = 0
    for node in root.get_nonterminals():
        if len(node.clades) != 2:
            raise ValueError("Expected binary tumour subtree")
        if rename or node.name is None:
            counter += 1
            while f"inferred_{counter}" in used:
                counter += 1
            node.name = f"inferred_{counter}"
            used.add(node.name)
    result = Tree(root=Clade(name=normal_name, branch_length=0, clades=[root]))
    names = [n.name for n in result.find_clades()]
    if len(names) != len(set(names)):
        raise ValueError("Node names must be unique")
    for node in result.find_clades():
        node.branch_length = 0
    return result


def infer_tree(distances, normal_name, out, executable="fastme"):
    """FastME BME, NNI and SPR with the seed used in the current experiments."""
    binary = shutil.which(str(executable))
    if binary is None:
        raise ValueError(f"FastME executable not found: {executable}")
    labels = list(distances.index)
    if any(re.search(r"[\s(),:;\[\]]", name) for name in labels):
        raise ValueError(
            "FastME sample IDs cannot contain whitespace or Newick punctuation"
        )
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    matrix = out / "distances.phy"
    width = max(map(len, labels)) + 1
    with matrix.open("w") as handle:
        handle.write(str(len(labels)) + "\n")
        for label, row in zip(labels, distances.to_numpy()):
            handle.write(label.ljust(width) + " ".join(f"{v:.5f}" for v in row) + "\n")
    command = [
        binary,
        "-i",
        str(matrix),
        "-o",
        str(out / "fastme_raw.nwk"),
        "-I",
        str(out / "fastme_stats.txt"),
        "-m",
        "B",
        "-z",
        "104000",
        "-T",
        "1",
        "-n",
        "-s",
    ]
    process = subprocess.run(command, cwd=out, text=True, capture_output=True)
    (out / "fastme.log").write_text(process.stdout + "\n" + process.stderr)
    process.check_returncode()
    tree = Phylo.read(out / "fastme_raw.nwk", "newick")
    if {c.name for c in tree.get_terminals()} != set(labels):
        raise ValueError("FastME leaf labels differ")
    return normalize_tree(tree, normal_name, labels, rename=True), command


def decode(profiles, coordinates, columns):
    sizes = []
    previous = None
    for chrom, _, _ in coordinates:
        if chrom != previous:
            sizes.append(0)
            previous = chrom
        sizes[-1] += 1
    frames = {}
    for name, word in profiles.items():
        chunks = word.split("X")
        if [len(chunk) for chunk in chunks] != sizes * len(columns):
            raise ValueError("Reconstructed bin layout changed")
        values = np.frombuffer("".join(chunks).encode(), dtype=np.uint8) - 48
        frames[name] = pd.DataFrame(
            values.reshape(len(columns), len(coordinates)).T,
            index=coordinates,
            columns=columns,
        )
    return pd.concat(frames, names=["sample_id"]).sort_index()


def reconstruct(
    input_df,
    normal_name="diploid",
    input_tree=None,
    n_cores=1,
    fastme="fastme",
    nni=False,
    topology_only=False,
    output_dir=None,
    allele_columns=None,
    nni_max_iter=math.inf,
    nni_max_frontier=1,
    nni_stop_on_plateau=True,
    progress=None,
):
    """Public reference-tool API. Input is the standard MEDICC copy-number dataframe.

    Returns a Result containing distances, rooted trees, profiles and diagnostics.
    No experiment manifest, truth profile or evaluation metric is needed.
    """
    if nni and topology_only:
        raise ValueError("NNI requires ancestor reconstruction")
    workers = n_cores or 1
    if workers < 1:
        raise ValueError("n_cores must be positive")
    genomes, coordinates, columns = encode(input_df, normal_name, allele_columns)
    labels = list(genomes)
    distances = None
    details = dict(model="full L_LOH", max_cn=8, max_wgds=2, prune_weight=0, nni=nni)
    if input_tree is None:
        distances, stats = distance_matrix(genomes, workers, progress)
        details["shortening_statistics"] = stats
        if output_dir is None:
            with tempfile.TemporaryDirectory(prefix="medicc-full-fst-") as tmp:
                tree, command = infer_tree(distances, normal_name, tmp, fastme)
        else:
            tree, command = infer_tree(distances, normal_name, output_dir, fastme)
        details["fastme_command"] = command
    else:
        if isinstance(input_tree, (str, Path)):
            input_tree = Phylo.read(input_tree, "newick")
        tree = normalize_tree(input_tree, normal_name, labels)
    initial = copy.deepcopy(tree)
    profiles = None
    search = None
    if nni:
        import medicc.core
        from .ancestors import FULL_TWO, Genome, update_lengths

        search = medicc.core.nni_mode(
            tree,
            {k: Genome(v) for k, v in genomes.items()},
            FULL_TWO,
            FULL_TWO,
            normal_name=normal_name,
            prune_weight=0,
            n_cores=workers,
            nni_max_iter=nni_max_iter,
            nni_max_frontier=nni_max_frontier,
            nni_stop_on_plateau=nni_stop_on_plateau,
            nni_trace_dir=str(output_dir) if output_dir is not None else None,
        )
        tree = search["best_trees"][0]
        profiles = search["best_ancestors"][0]
        if profiles is None:
            from .ancestors import reconstruct as resolve_tree

            profiles = resolve_tree(
                tree,
                {k: Genome(v) for k, v in genomes.items()},
                FULL_TWO,
                FULL_TWO,
                normal_name,
                n_cores=workers,
            )
        update_lengths(tree, FULL_TWO, profiles, normal_name)
        if (
            sum(n.branch_length or 0 for n in tree.find_clades())
            != search["best_score"]
        ):
            raise ValueError("Final directed scores disagree with NNI")
        details.update(
            initial_score=search["trace"][0],
            final_score=search["best_score"],
            score_trace=search["trace"],
        )
    elif not topology_only:
        root = tree.root.clades[0]
        children = {n.name: [c.name for c in n.clades] for n in root.get_nonterminals()}
        leaves = {k: v for k, v in genomes.items() if k != normal_name}
        profiles, _, metrics = reconstruct_messages(
            children,
            leaves,
            root.name,
            genomes[normal_name],
            workers=workers,
            state_limit=0,
            progress=progress,
        )
        profiles[normal_name] = genomes[normal_name]
        model = Model()
        try:
            for parent in tree.find_clades():
                for child in parent.clades:
                    value = model.scorer.score(
                        *model.shorten(profiles[parent.name], profiles[child.name]),
                        directed=True,
                    )
                    if not math.isfinite(value):
                        raise ValueError("Infeasible reconstructed branch")
                    child.branch_length = value
        finally:
            model.close()
        details.update(metrics)
    output_df = (
        decode(profiles, coordinates, columns)
        if profiles is not None
        else input_df.copy()
    )
    details["tree_length"] = sum(n.branch_length or 0 for n in tree.find_clades())
    return Result(
        labels, distances, initial, tree, profiles, output_df, details, search
    )
