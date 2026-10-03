#!/usr/bin/env python3

"""
Measure structural diversity of RFD3-generated cyclic macrocycles.

For each peptide length, the script:

1. Finds all generated .cif.gz structures.
2. Creates one deterministic shuffled sampling order.
3. Saves that order to sampling_order.csv.
4. Clusters nested prefixes of that order.

By default, the sampled fractions are:

    1/8, 1/4, 1/2, 1

For 10,000 structures this gives:

    1,250
    2,500
    5,000
    10,000

This mirrors the sampling scheme used in RFpeptides Supplementary Fig. S2H,
where 48,000 structures were sampled at approximately:

    6,000
    12,000
    24,000
    48,000

Clustering uses PyRosetta EnergyBasedClusteringProtocol with:

    - backbone Cartesian RMSD
    - 0.5 Angstrom clustering radius
    - cyclic peptide handling
    - all cyclic residue permutations
    - no C-beta atoms
    - no prerelaxation

Important:
Rosetta's bb_cartesian metric is backbone Cartesian RMSD rather than strictly
C-alpha-only RMSD. This is therefore analogous to, but not an atom-for-atom
reproduction of, the RFpeptides clustering metric reported for figure S2H.

Outputs for each peptide length:

    length_10/
        sampling_order.csv
        cluster_growth.csv
        clustering_result.csv
        summary.json
        cluster_centers/
            manifest.csv
            summary.json
            *.cif.gz

    length_12/
        sampling_order.csv
        cluster_growth.csv
        clustering_result.csv
        summary.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path


# =============================================================================
# Paths
# =============================================================================

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1]

DEFAULT_INPUT_DIR = (
    EXPERIMENT_DIR
    / "test_2"
    / "macrocycle_monomer_10k"
)

DEFAULT_OUTPUT_DIR = (
    DEFAULT_INPUT_DIR
    / "clustering"
)


# =============================================================================
# Filename parsing
# =============================================================================

# Example:
#
# macrocycle_monomer_10K_macrocycle_10_7_model_314.cif.gz
#
#                                     ^^ batch
#                                             ^^^ model
#
INPUT_FILENAME_PATTERN = re.compile(
    r"macrocycle_(?P<length>\d+)"
    r"_(?P<batch>\d+)"
    r"_model_(?P<model>\d+)"
    r"\.cif\.gz$"
)

CLUSTER_CENTER_LOG_PATTERN = re.compile(
    r"Started cluster (?P<cluster_id>\d+) "
    r"and added structure (?P<structure_index>\d+) to it\."
)


# =============================================================================
# Data structures
# =============================================================================


@dataclass(frozen=True)
class DesignFile:
    """One generated macrocycle structure."""

    path: Path
    length: int
    batch: int
    model: int

    @property
    def generation_order(self) -> tuple[int, int]:
        return self.batch, self.model


@dataclass(frozen=True)
class ClusteringResult:
    """Result from clustering one sampled subset."""

    peptide_length: int
    n_structures: int
    n_clusters: int
    unique_cluster_fraction: float
    elapsed_seconds: float


# =============================================================================
# Arguments
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cluster cyclic macrocycle backbones with PyRosetta."
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing .cif.gz structures. Default: {DEFAULT_INPUT_DIR}",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for clustering results. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=[10, 12],
        help="Peptide lengths to analyze. Default: 10 12",
    )

    parser.add_argument(
        "--radius",
        type=float,
        default=0.5,
        help="Backbone Cartesian RMSD clustering radius in Angstroms. Default: 0.5",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help=(
            "Base random seed used to construct deterministic sampling orders. "
            "Default: 2026"
        ),
    )

    parser.add_argument(
        "--sample-sizes",
        type=int,
        nargs="+",
        default=None,
        help=(
            "Optional explicit sample sizes. "
            "Default: 1/8, 1/4, 1/2, and all available structures."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing results for the requested peptide lengths.",
    )

    parser.add_argument(
        "--save-cluster-centers-only",
        action="store_true",
        help=(
            "Reuse each existing sampling_order.csv, rerun only the full "
            "dataset calculation, and save one original input structure per "
            "cluster. Existing cluster-growth results remain unchanged."
        ),
    )

    # -------------------------------------------------------------------------
    # Internal worker arguments.
    #
    # Each clustering calculation gets a fresh process because PyRosetta is
    # most safely initialized once per process.
    # -------------------------------------------------------------------------

    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_worker-length", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--_sample-size", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--_sampling-order", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--_result-path", type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--_save-cluster-centers",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--_center-output-dir",
        type=Path,
        help=argparse.SUPPRESS,
    )

    return parser.parse_args()


# =============================================================================
# Input discovery
# =============================================================================


def discover_designs(
    input_dir: Path,
    peptide_length: int,
) -> list[DesignFile]:
    """Find all structures for one peptide length."""

    designs: list[DesignFile] = []

    for path in input_dir.glob("*.cif.gz"):
        match = INPUT_FILENAME_PATTERN.search(path.name)

        if match is None:
            continue

        design = DesignFile(
            path=path.resolve(),
            length=int(match.group("length")),
            batch=int(match.group("batch")),
            model=int(match.group("model")),
        )

        if design.length == peptide_length:
            designs.append(design)

    designs.sort(key=lambda design: design.generation_order)

    return designs


def validate_designs(designs: list[DesignFile]) -> None:
    """Check that batch/model identifiers are unique."""

    seen: set[tuple[int, int]] = set()

    for design in designs:
        identifier = (design.batch, design.model)

        if identifier in seen:
            raise ValueError(
                "Duplicate structure identifier: "
                f"batch={design.batch}, model={design.model}"
            )

        seen.add(identifier)


# =============================================================================
# Sampling
# =============================================================================


def determine_sample_sizes(
    total_structures: int,
    requested_sizes: list[int] | None,
) -> list[int]:
    """
    Determine which nested sample sizes will be clustered.

    Default:
        1/8, 1/4, 1/2, and all structures.

    The complete dataset is always included.
    """

    if requested_sizes is None:
        sample_sizes = [
            max(1, total_structures // 8),
            max(1, total_structures // 4),
            max(1, total_structures // 2),
            total_structures,
        ]
    else:
        sample_sizes = list(requested_sizes)

        for sample_size in sample_sizes:
            if sample_size <= 0:
                raise ValueError("Sample sizes must be greater than zero.")

            if sample_size > total_structures:
                raise ValueError(
                    f"Requested sample size {sample_size:,} exceeds the "
                    f"{total_structures:,} available structures."
                )

        # The final/full result is always useful and remains our canonical
        # clustering result.
        sample_sizes.append(total_structures)

    return sorted(set(sample_sizes))


def make_sampling_order(
    designs: list[DesignFile],
    base_seed: int,
    peptide_length: int,
) -> tuple[list[DesignFile], int]:
    """
    Shuffle structures once and deterministically.

    A different deterministic seed is used for each peptide length so the
    10-mer and 12-mer populations are sampled independently.

    All sample sizes are prefixes of this same shuffled list.
    """

    effective_seed = base_seed + peptide_length

    shuffled = list(designs)

    rng = random.Random(effective_seed)
    rng.shuffle(shuffled)

    return shuffled, effective_seed


def write_sampling_order(
    designs: list[DesignFile],
    output_path: Path,
) -> None:
    """Save the exact shuffled order used for nested sampling."""

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow(
            [
                "sample_rank",
                "peptide_length",
                "batch",
                "model",
                "filename",
                "original_path",
            ]
        )

        for rank, design in enumerate(designs, start=1):
            writer.writerow(
                [
                    rank,
                    design.length,
                    design.batch,
                    design.model,
                    design.path.name,
                    str(design.path),
                ]
            )


def read_sampling_order(
    sampling_order_path: Path,
    sample_size: int,
) -> list[DesignFile]:
    """Read the first N structures from a saved sampling order."""

    designs: list[DesignFile] = []

    with sampling_order_path.open() as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            designs.append(
                DesignFile(
                    path=Path(row["original_path"]),
                    length=int(row["peptide_length"]),
                    batch=int(row["batch"]),
                    model=int(row["model"]),
                )
            )

            if len(designs) == sample_size:
                break

    if len(designs) != sample_size:
        raise ValueError(
            f"Requested {sample_size:,} structures but only "
            f"{len(designs):,} were available in {sampling_order_path}."
        )

    return designs


# =============================================================================
# Rosetta input preparation
# =============================================================================


def decompress_designs(
    designs: list[DesignFile],
    scratch_dir: Path,
) -> list[Path]:
    """Decompress the selected .cif.gz structures into a temporary directory."""

    cif_paths: list[Path] = []
    total = len(designs)

    print(f"Decompressing {total:,} mmCIF files...", flush=True)

    for index, design in enumerate(designs, start=1):
        output_path = scratch_dir / design.path.name.removesuffix(".gz")

        with gzip.open(design.path, "rb") as source:
            with output_path.open("wb") as destination:
                shutil.copyfileobj(source, destination)

        cif_paths.append(output_path)

        if index % 1000 == 0 or index == total:
            print(f"  {index:,}/{total:,}", flush=True)

    return cif_paths


def write_rosetta_input_files(
    cif_paths: list[Path],
    scratch_dir: Path,
) -> tuple[Path, Path]:
    """
    Write Rosetta's structure list and deterministic external score file.

    Scores are simply the position in the shuffled sampling order:

        first structure  -> 0
        second structure -> 1
        third structure  -> 2

    Because every larger sample is a prefix of the same ordering, existing
    cluster centers retain the same priority as sample size increases.
    """

    input_list_path = scratch_dir / "rosetta_inputs.txt"
    score_file_path = scratch_dir / "alternative_scores.txt"

    with input_list_path.open("w") as input_handle:
        with score_file_path.open("w") as score_handle:
            for rank, cif_path in enumerate(cif_paths):
                input_handle.write(f"{cif_path.name}\n")
                score_handle.write(f"{cif_path.name} {float(rank)}\n")

    return input_list_path, score_file_path


def read_cluster_center_indices(
    log_path: Path,
    *,
    expected_cluster_count: int,
    structure_count: int,
) -> list[tuple[int, int]]:
    """Read ``(cluster ID, structure index)`` pairs from Rosetta's log."""

    centers: list[tuple[int, int]] = []

    with log_path.open() as handle:
        for line in handle:
            match = CLUSTER_CENTER_LOG_PATTERN.search(line)

            if match is None:
                continue

            centers.append(
                (
                    int(match.group("cluster_id")),
                    int(match.group("structure_index")),
                )
            )

    if len(centers) != expected_cluster_count:
        raise RuntimeError(
            f"Rosetta reported {expected_cluster_count:,} clusters but "
            f"{len(centers):,} center records were found in {log_path}."
        )

    expected_cluster_ids = list(range(1, expected_cluster_count + 1))
    observed_cluster_ids = [cluster_id for cluster_id, _ in centers]

    if observed_cluster_ids != expected_cluster_ids:
        raise RuntimeError(
            "Cluster IDs in the Rosetta log are missing or out of order."
        )

    for _, structure_index in centers:
        if not 1 <= structure_index <= structure_count:
            raise RuntimeError(
                f"Rosetta center index {structure_index} is outside the "
                f"1..{structure_count} input range."
            )

    return centers


def save_cluster_center_structures(
    *,
    centers: list[tuple[int, int]],
    designs: list[DesignFile],
    output_dir: Path,
) -> Path:
    """Copy the original generated structure selected for each cluster."""

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "manifest.csv"

    with manifest_path.open("w", newline="") as handle:
        fieldnames = [
            "cluster_id",
            "sampling_rank",
            "peptide_length",
            "batch",
            "model",
            "filename",
            "source_path",
            "representative_path",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()

        for cluster_id, structure_index in centers:
            design = designs[structure_index - 1]
            representative_path = output_dir / design.path.name

            shutil.copy2(design.path, representative_path)

            writer.writerow(
                {
                    "cluster_id": cluster_id,
                    "sampling_rank": structure_index,
                    "peptide_length": design.length,
                    "batch": design.batch,
                    "model": design.model,
                    "filename": design.path.name,
                    "source_path": str(design.path),
                    "representative_path": str(
                        representative_path.resolve()
                    ),
                }
            )

    return manifest_path


# =============================================================================
# PyRosetta clustering
# =============================================================================


def run_rosetta_clustering(
    input_list_path: Path,
    score_file_path: Path,
    scratch_dir: Path,
    radius: float,
    save_cluster_centers: bool,
) -> tuple[int, str, float, list[tuple[int, int]]]:
    """Run one complete Rosetta clustering calculation."""

    import pyrosetta
    from pyrosetta.rosetta.protocols import energy_based_clustering as ebc

    original_working_directory = Path.cwd()

    try:
        os.chdir(scratch_dir)

        initialization_options = [
            f"-in:file:l {input_list_path.name}",
            "-in:file:fullatom",
            "-mute all",
        ]

        tracer_log_path = scratch_dir / "rosetta_clustering.log"

        if save_cluster_centers:
            tracer_channel = (
                "protocols.cluster.energy_based_clustering."
                "EnergyBasedClusteringProtocol"
            )
            initialization_options.extend(
                [
                    "-unmute",
                    tracer_channel,
                ]
            )

            rosetta_logger = logging.getLogger("rosetta")
            rosetta_logger.setLevel(logging.INFO)
            rosetta_logger.propagate = False
            rosetta_logger.addHandler(
                logging.FileHandler(
                    tracer_log_path,
                    mode="w",
                )
            )

            pyrosetta.init(
                " ".join(initialization_options),
                set_logging_handler="logging",
            )
        else:
            pyrosetta.init(" ".join(initialization_options))

        options = ebc.EnergyBasedClusteringOptions(False)

        # Structural distance metric.
        options.cluster_by_ = ebc.EBC_bb_cartesian
        options.cluster_radius_ = float(radius)
        options.use_CB_ = False

        # Cyclic-peptide comparison.
        options.cyclic_ = True
        options.cluster_cyclic_permutations_ = True
        options.cyclic_permutation_offset_ = 1

        # Preserve the generated structures exactly.
        options.prerelax_ = False
        options.mutate_to_ala_ = False

        # Use sampling order, rather than Rosetta energy, to choose cluster
        # centers.
        options.path_to_scores_file_ = str(score_file_path.resolve())

        # Rosetta still performs the complete clustering calculation
        # internally. Usually only one diagnostic output is requested. When
        # saving centers, one output per cluster is needed so every center is
        # present in the tracer log and Rosetta's own output table.
        options.limit_structures_per_cluster_ = 1
        options.limit_clusters_ = 0 if save_cluster_centers else 1
        options.silent_output_ = False

        # Any diagnostic PDB written by Rosetta stays inside the temporary
        # scratch directory and disappears after the run.
        options.output_prefix_ = str(scratch_dir.resolve() / "rosetta_")

        protocol = ebc.EnergyBasedClusteringProtocol(options)

        print("Running Rosetta clustering...", flush=True)

        start_time = time.time()
        protocol.go()
        elapsed_seconds = time.time() - start_time

        n_clusters = int(protocol.n_clusters_from_last_run())

        if n_clusters <= 0:
            raise RuntimeError("Rosetta completed but reported zero clusters.")

        try:
            pyrosetta_version = str(pyrosetta.version())
        except Exception:
            pyrosetta_version = "unknown"

        print(
            f"Rosetta finished: {n_clusters:,} clusters "
            f"in {elapsed_seconds:.2f} seconds.",
            flush=True,
        )

        center_indices: list[tuple[int, int]] = []

        if save_cluster_centers:
            for handler in logging.getLogger("rosetta").handlers:
                handler.flush()

            center_indices = read_cluster_center_indices(
                tracer_log_path,
                expected_cluster_count=n_clusters,
                structure_count=sum(
                    1
                    for line in input_list_path.read_text().splitlines()
                    if line.strip()
                ),
            )

        return (
            n_clusters,
            pyrosetta_version,
            elapsed_seconds,
            center_indices,
        )

    finally:
        os.chdir(original_working_directory)


# =============================================================================
# Worker process
# =============================================================================


def worker_main(args: argparse.Namespace) -> None:
    """
    Cluster one sample size.

    Each worker initializes PyRosetta exactly once.
    """

    required = {
        "--_worker-length": args._worker_length,
        "--_sample-size": args._sample_size,
        "--_sampling-order": args._sampling_order,
        "--_result-path": args._result_path,
    }

    missing = [name for name, value in required.items() if value is None]

    if args._save_cluster_centers and args._center_output_dir is None:
        missing.append("--_center-output-dir")

    if missing:
        raise ValueError(
            "Worker mode is missing required arguments: "
            + ", ".join(missing)
        )

    peptide_length = args._worker_length
    sample_size = args._sample_size
    sampling_order_path = args._sampling_order.resolve()
    result_path = args._result_path.resolve()

    designs = read_sampling_order(
        sampling_order_path,
        sample_size,
    )

    scratch_root = args.output_dir.resolve() / "_scratch"
    scratch_root.mkdir(parents=True, exist_ok=True)

    print()
    print("-" * 72, flush=True)
    print(
        f"{peptide_length}-mer | sampled structures: {sample_size:,}",
        flush=True,
    )
    print("-" * 72, flush=True)

    with tempfile.TemporaryDirectory(
        prefix=f"length_{peptide_length}_n{sample_size}_",
        dir=scratch_root,
    ) as temporary_directory:
        scratch_dir = Path(temporary_directory)

        cif_paths = decompress_designs(
            designs,
            scratch_dir,
        )

        input_list_path, score_file_path = write_rosetta_input_files(
            cif_paths,
            scratch_dir,
        )

        (
            n_clusters,
            pyrosetta_version,
            elapsed_seconds,
            center_indices,
        ) = run_rosetta_clustering(
            input_list_path=input_list_path,
            score_file_path=score_file_path,
            scratch_dir=scratch_dir,
            radius=args.radius,
            save_cluster_centers=args._save_cluster_centers,
        )

        manifest_path = None

        if args._save_cluster_centers:
            manifest_path = save_cluster_center_structures(
                centers=center_indices,
                designs=designs,
                output_dir=args._center_output_dir.resolve(),
            )

    result = ClusteringResult(
        peptide_length=peptide_length,
        n_structures=sample_size,
        n_clusters=n_clusters,
        unique_cluster_fraction=n_clusters / sample_size,
        elapsed_seconds=elapsed_seconds,
    )

    result_path.parent.mkdir(parents=True, exist_ok=True)

    with result_path.open("w") as handle:
        json.dump(
            {
                **asdict(result),
                "pyrosetta_version": pyrosetta_version,
                "cluster_center_count": len(center_indices),
                "cluster_center_manifest": (
                    str(manifest_path)
                    if manifest_path is not None
                    else None
                ),
            },
            handle,
            indent=2,
        )
        handle.write("\n")

    try:
        scratch_root.rmdir()
    except OSError:
        pass


# =============================================================================
# Permanent output
# =============================================================================


def write_cluster_growth(
    results: list[ClusteringResult],
    output_path: Path,
) -> None:
    """Write all sampled points for the growth curve."""

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "peptide_length",
                "n_structures",
                "n_clusters",
                "unique_cluster_fraction",
                "elapsed_seconds",
            ],
        )

        writer.writeheader()

        for result in results:
            writer.writerow(asdict(result))


def write_final_result(
    result: ClusteringResult,
    output_path: Path,
) -> None:
    """
    Preserve a simple one-row final result for convenience.

    This is the full-dataset point only.
    """

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "peptide_length",
                "n_structures",
                "n_clusters",
                "unique_cluster_fraction",
                "elapsed_seconds",
            ],
        )

        writer.writeheader()
        writer.writerow(asdict(result))


def write_summary(
    *,
    output_path: Path,
    peptide_length: int,
    total_structures: int,
    sample_sizes: list[int],
    base_seed: int,
    effective_seed: int,
    radius: float,
    pyrosetta_version: str,
    results: list[ClusteringResult],
    input_dir: Path,
) -> None:
    """Write method settings and run provenance."""

    summary = {
        "peptide_length": peptide_length,
        "total_structures": total_structures,
        "sample_sizes": sample_sizes,
        "sampling": {
            "method": "nested prefixes of one deterministic shuffled ordering",
            "base_seed": base_seed,
            "effective_seed": effective_seed,
        },
        "input_directory": str(input_dir.resolve()),
        "method": {
            "software": "PyRosetta EnergyBasedClusteringProtocol",
            "pyrosetta_version": pyrosetta_version,
            "cluster_by": "bb_cartesian",
            "cluster_radius_angstrom": radius,
            "use_CB": False,
            "cyclic": True,
            "cluster_cyclic_permutations": True,
            "cyclic_permutation_offset": 1,
            "prerelax": False,
            "mutate_to_ala": False,
            "seed_order": (
                "sampling rank supplied as deterministic external score"
            ),
            "rosetta_output": {
                "limit_structures_per_cluster": 1,
                "limit_clusters": 1,
                "note": (
                    "These options limit files written to disk, not the "
                    "complete clustering calculation."
                ),
            },
        },
        "results": [asdict(result) for result in results],
    }

    with output_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")


# =============================================================================
# Directory handling
# =============================================================================


def prepare_length_directory(
    output_dir: Path,
    peptide_length: int,
    overwrite: bool,
) -> Path:
    """Create a clean permanent directory for one peptide length."""

    length_dir = output_dir / f"length_{peptide_length}"

    if length_dir.exists():
        if not overwrite and any(length_dir.iterdir()):
            raise FileExistsError(
                f"Output already exists:\n"
                f"  {length_dir}\n\n"
                "Use --overwrite to replace it."
            )

        if overwrite:
            shutil.rmtree(length_dir)

    length_dir.mkdir(parents=True, exist_ok=True)

    return length_dir


# =============================================================================
# Parent process
# =============================================================================


def analyze_length(
    args: argparse.Namespace,
    peptide_length: int,
) -> None:
    """Run the complete cluster-growth analysis for one peptide length."""

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    designs = discover_designs(
        input_dir,
        peptide_length,
    )

    if not designs:
        raise FileNotFoundError(
            f"No {peptide_length}-mer .cif.gz files found in:\n{input_dir}"
        )

    validate_designs(designs)

    total_structures = len(designs)

    sample_sizes = determine_sample_sizes(
        total_structures,
        args.sample_sizes,
    )

    sampling_order, effective_seed = make_sampling_order(
        designs,
        base_seed=args.seed,
        peptide_length=peptide_length,
    )

    length_dir = prepare_length_directory(
        output_dir,
        peptide_length,
        args.overwrite,
    )

    sampling_order_path = length_dir / "sampling_order.csv"

    write_sampling_order(
        sampling_order,
        sampling_order_path,
    )

    print()
    print("=" * 72, flush=True)
    print(f"Analyzing {peptide_length}-mer macrocycles", flush=True)
    print("=" * 72, flush=True)
    print(f"Available structures : {total_structures:,}", flush=True)
    print(
        "Sample sizes         : "
        + ", ".join(f"{value:,}" for value in sample_sizes),
        flush=True,
    )
    print(f"Sampling seed        : {effective_seed}", flush=True)
    print(f"Cluster radius       : {args.radius:.3f} Å", flush=True)
    print(f"Output directory     : {length_dir}", flush=True)

    results: list[ClusteringResult] = []
    pyrosetta_version = "unknown"

    # Worker results are temporary. Only the final clean CSV/JSON products
    # survive.
    with tempfile.TemporaryDirectory(
        prefix=f"_worker_results_length_{peptide_length}_",
        dir=output_dir,
    ) as temporary_directory:
        worker_result_dir = Path(temporary_directory)

        for sample_size in sample_sizes:
            result_path = (
                worker_result_dir
                / f"n{sample_size}.json"
            )

            command = [
                sys.executable,
                str(SCRIPT_PATH),
                "--_worker",
                "--_worker-length",
                str(peptide_length),
                "--_sample-size",
                str(sample_size),
                "--_sampling-order",
                str(sampling_order_path),
                "--_result-path",
                str(result_path),
                "--output-dir",
                str(output_dir),
                "--radius",
                str(args.radius),
            ]

            subprocess.run(
                command,
                check=True,
            )

            with result_path.open() as handle:
                worker_output = json.load(handle)

            pyrosetta_version = worker_output["pyrosetta_version"]

            results.append(
                ClusteringResult(
                    peptide_length=int(worker_output["peptide_length"]),
                    n_structures=int(worker_output["n_structures"]),
                    n_clusters=int(worker_output["n_clusters"]),
                    unique_cluster_fraction=float(
                        worker_output["unique_cluster_fraction"]
                    ),
                    elapsed_seconds=float(worker_output["elapsed_seconds"]),
                )
            )

    results.sort(key=lambda result: result.n_structures)

    # Sanity check: with nested prefixes and fixed cluster-center ordering,
    # discovered cluster count should never decrease.
    for previous, current in zip(results, results[1:]):
        if current.n_clusters < previous.n_clusters:
            raise RuntimeError(
                "Cluster count decreased as sample size increased:\n"
                f"  {previous.n_structures:,} -> {previous.n_clusters:,}\n"
                f"  {current.n_structures:,} -> {current.n_clusters:,}\n"
                "This should not happen with the nested deterministic sampling "
                "scheme and should be investigated."
            )

    write_cluster_growth(
        results,
        length_dir / "cluster_growth.csv",
    )

    final_result = results[-1]

    write_final_result(
        final_result,
        length_dir / "clustering_result.csv",
    )

    write_summary(
        output_path=length_dir / "summary.json",
        peptide_length=peptide_length,
        total_structures=total_structures,
        sample_sizes=sample_sizes,
        base_seed=args.seed,
        effective_seed=effective_seed,
        radius=args.radius,
        pyrosetta_version=pyrosetta_version,
        results=results,
        input_dir=input_dir,
    )

    print()
    print(f"{peptide_length}-mer cluster growth:", flush=True)

    for result in results:
        print(
            f"  {result.n_structures:>10,} structures -> "
            f"{result.n_clusters:>6,} clusters "
            f"({result.unique_cluster_fraction:.2%})",
            flush=True,
        )


def parent_main(args: argparse.Namespace) -> None:
    """Run the requested peptide lengths sequentially."""

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist:\n{input_dir}"
        )

    if args.radius <= 0:
        raise ValueError("--radius must be greater than zero.")

    output_dir.mkdir(parents=True, exist_ok=True)

    for peptide_length in args.lengths:
        analyze_length(
            args,
            peptide_length,
        )

    print()
    print("=" * 72, flush=True)
    print("Clustering analysis complete.", flush=True)
    print(f"Results: {output_dir}", flush=True)
    print("=" * 72, flush=True)


# =============================================================================
# Cluster-center-only mode
# =============================================================================


def read_canonical_cluster_count(
    result_path: Path,
    expected_structure_count: int,
) -> int:
    """Read the existing full-dataset cluster count."""

    with result_path.open() as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 1:
        raise ValueError(
            f"Expected one row in {result_path}; found {len(rows)}."
        )

    structure_count = int(rows[0]["n_structures"])

    if structure_count != expected_structure_count:
        raise ValueError(
            f"Existing result uses {structure_count:,} structures, but "
            f"{expected_structure_count:,} are currently available."
        )

    return int(rows[0]["n_clusters"])


def prepare_cluster_center_directory(
    center_dir: Path,
    *,
    overwrite: bool,
) -> None:
    """Create an empty directory for permanent cluster representatives."""

    if center_dir.exists() and any(center_dir.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"Cluster-center output already exists:\n  {center_dir}\n\n"
                "Use --overwrite to replace it."
            )

        shutil.rmtree(center_dir)

    center_dir.mkdir(parents=True, exist_ok=True)


def save_cluster_centers_for_length(
    args: argparse.Namespace,
    peptide_length: int,
) -> None:
    """Rerun one full clustering point and save its representatives."""

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    length_dir = output_dir / f"length_{peptide_length}"
    sampling_order_path = length_dir / "sampling_order.csv"
    canonical_result_path = length_dir / "clustering_result.csv"
    canonical_summary_path = length_dir / "summary.json"

    for required_path in (
        sampling_order_path,
        canonical_result_path,
        canonical_summary_path,
    ):
        if not required_path.is_file():
            raise FileNotFoundError(
                "Cluster-center-only mode requires the existing clustering "
                f"output:\n{required_path}"
            )

    designs = discover_designs(input_dir, peptide_length)
    validate_designs(designs)
    total_structures = len(designs)

    if total_structures == 0:
        raise FileNotFoundError(
            f"No {peptide_length}-mer structures found in:\n{input_dir}"
        )

    # Reading the complete saved order verifies that it still references the
    # same number of available structures before a costly clustering run.
    ordered_designs = read_sampling_order(
        sampling_order_path,
        total_structures,
    )

    missing_sources = [
        design.path
        for design in ordered_designs
        if not design.path.is_file()
    ]

    if missing_sources:
        raise FileNotFoundError(
            "A structure referenced by sampling_order.csv is missing:\n"
            f"{missing_sources[0]}"
        )

    expected_cluster_count = read_canonical_cluster_count(
        canonical_result_path,
        total_structures,
    )

    with canonical_summary_path.open() as handle:
        canonical_summary = json.load(handle)

    canonical_radius = float(
        canonical_summary["method"]["cluster_radius_angstrom"]
    )

    if args.radius != canonical_radius:
        raise ValueError(
            "The requested radius does not match the existing clustering "
            f"result: requested {args.radius}, existing {canonical_radius}."
        )

    center_dir = length_dir / "cluster_centers"
    prepare_cluster_center_directory(
        center_dir,
        overwrite=args.overwrite,
    )

    print()
    print("=" * 72, flush=True)
    print(
        f"Saving full-dataset {peptide_length}-mer cluster centers",
        flush=True,
    )
    print("=" * 72, flush=True)
    print(f"Structures             : {total_structures:,}", flush=True)
    print(f"Expected clusters      : {expected_cluster_count:,}", flush=True)
    print(f"Saved sampling order   : {sampling_order_path}", flush=True)
    print(f"Center output directory: {center_dir}", flush=True)

    with tempfile.TemporaryDirectory(
        prefix=f"_cluster_centers_length_{peptide_length}_",
        dir=output_dir,
    ) as temporary_directory:
        result_path = Path(temporary_directory) / "result.json"

        command = [
            sys.executable,
            str(SCRIPT_PATH),
            "--_worker",
            "--_worker-length",
            str(peptide_length),
            "--_sample-size",
            str(total_structures),
            "--_sampling-order",
            str(sampling_order_path),
            "--_result-path",
            str(result_path),
            "--_save-cluster-centers",
            "--_center-output-dir",
            str(center_dir),
            "--output-dir",
            str(output_dir),
            "--radius",
            str(args.radius),
        ]

        subprocess.run(command, check=True)

        with result_path.open() as handle:
            worker_output = json.load(handle)

    observed_cluster_count = int(worker_output["n_clusters"])
    saved_center_count = int(worker_output["cluster_center_count"])

    if observed_cluster_count != expected_cluster_count:
        raise RuntimeError(
            "The repeated full-dataset clustering did not reproduce the "
            "canonical cluster count:\n"
            f"  expected: {expected_cluster_count:,}\n"
            f"  observed: {observed_cluster_count:,}"
        )

    if saved_center_count != expected_cluster_count:
        raise RuntimeError(
            f"Expected {expected_cluster_count:,} saved centers but found "
            f"{saved_center_count:,}."
        )

    summary = {
        "peptide_length": peptide_length,
        "structure_count": total_structures,
        "cluster_count": observed_cluster_count,
        "cluster_radius_angstrom": args.radius,
        "elapsed_seconds": float(worker_output["elapsed_seconds"]),
        "pyrosetta_version": worker_output["pyrosetta_version"],
        "sampling_order": str(sampling_order_path),
        "manifest": worker_output["cluster_center_manifest"],
        "representative_coordinates": (
            "unchanged copies of the original generated mmCIF structures "
            "selected as cluster centers"
        ),
        "center_selection": (
            "Rosetta structure index parsed from the deterministic "
            "EnergyBasedClusteringProtocol tracer log"
        ),
    }

    summary_path = center_dir / "summary.json"

    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    ignore_path = center_dir / ".gitignore"
    ignore_path.write_text(
        "# Reproducible copies listed in manifest.csv\n"
        "*.cif.gz\n"
    )

    print(
        f"Saved {saved_center_count:,} representatives to {center_dir}",
        flush=True,
    )


def cluster_centers_only_main(args: argparse.Namespace) -> None:
    """Save representatives without replacing the cluster-growth outputs."""

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist:\n{input_dir}"
        )

    if not output_dir.is_dir():
        raise NotADirectoryError(
            "Existing clustering directory does not exist:\n"
            f"{output_dir}"
        )

    if args.radius <= 0:
        raise ValueError("--radius must be greater than zero.")

    for peptide_length in args.lengths:
        save_cluster_centers_for_length(args, peptide_length)

    print()
    print("Cluster-center export complete.", flush=True)
    print(f"Results: {output_dir}", flush=True)


# =============================================================================
# Entry point
# =============================================================================


def main() -> None:
    args = parse_arguments()

    if args._worker:
        worker_main(args)
    elif args.save_cluster_centers_only:
        cluster_centers_only_main(args)
    else:
        parent_main(args)


if __name__ == "__main__":
    main()
