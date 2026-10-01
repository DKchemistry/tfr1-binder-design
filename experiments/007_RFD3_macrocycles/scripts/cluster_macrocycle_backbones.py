#!/usr/bin/env python3

"""
Cluster RFD3-generated cyclic macrocycle backbones with PyRosetta.

This script is intended to measure backbone diversity in a way analogous to
the RFpeptides Supplementary Fig. S2H analysis.

For each peptide length (10 and 12 residues by default), it:

1. Finds all matching RFD3 .cif.gz structures.
2. Sorts them deterministically by batch number and model number.
3. Temporarily decompresses them to ordinary .cif files.
4. Supplies external scores equal to input rank so Rosetta uses a deterministic
   cluster-seed order rather than Rosetta energies.
5. Runs Rosetta EnergyBasedClusteringProtocol with:
       - backbone Cartesian RMSD
       - 0.5 Angstrom cluster radius
       - cyclic peptide handling
       - all one-residue cyclic permutations
       - no C-beta atoms
       - no prerelaxation
6. Records the total number of clusters returned by Rosetta.
7. Writes a CSV and JSON summary for downstream plotting.

Important
---------
Rosetta's "bb_cartesian" metric uses backbone Cartesian coordinates. It is
therefore closely analogous to, but not exactly identical to, the C-alpha-only
0.5 Angstrom clustering described in RFpeptides Supplementary Fig. S2H.

Why Rosetta output is deliberately limited
-------------------------------------------
EnergyBasedClusteringProtocol can reconstruct and superimpose cyclically
permuted cluster members when writing them to disk. In the PyRosetta build
used for this analysis, that output path can trigger a numerical assertion
inside Rosetta's rms_util.cc.

The scientific quantity needed here is the TOTAL NUMBER OF CLUSTERS, which
Rosetta computes before applying its output limits. Therefore this script
asks Rosetta to write only one structure from one cluster. This avoids
unnecessary cyclic-member output alignment while preserving the complete
clustering calculation.

10-mers and 12-mers are run in separate Python subprocesses because Rosetta's
global options are most safely initialized once per process.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


# =============================================================================
# Default paths
# =============================================================================

SCRIPT_PATH = Path(__file__).resolve()

# Expected location:
#
# experiments/
#   007_RFD3_macrocycles/
#     scripts/
#       cluster_macrocycle_backbones.py
#
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
# Input filename parsing
# =============================================================================

# Example:
#
# macrocycle_monomer_10K_macrocycle_10_7_model_314.cif.gz
#                                      ^      ^^^
#                                    batch    model
#
INPUT_FILENAME_PATTERN = re.compile(
    r"macrocycle_(?P<length>\d+)"
    r"_(?P<batch>\d+)"
    r"_model_(?P<model>\d+)"
    r"\.cif\.gz$"
)


@dataclass(frozen=True)
class DesignFile:
    """Information parsed from one generated macrocycle structure."""

    path: Path
    length: int
    batch: int
    model: int

    @property
    def sort_key(self) -> tuple[int, int]:
        """Generation order: first batch number, then model number."""
        return self.batch, self.model


# =============================================================================
# Command-line arguments
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cluster cyclic RFD3 macrocycle backbones with PyRosetta."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=(
            "Directory containing RFD3 .cif.gz files. "
            f"Default: {DEFAULT_INPUT_DIR}"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Root directory for clustering results. "
            f"Default: {DEFAULT_OUTPUT_DIR}"
        ),
    )

    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=[10, 12],
        help="Peptide lengths to cluster separately. Default: 10 12",
    )

    parser.add_argument(
        "--radius",
        type=float,
        default=0.5,
        help=(
            "Backbone Cartesian RMSD clustering radius in Angstroms. "
            "Default: 0.5"
        ),
    )

    parser.add_argument(
        "--max-structures",
        type=int,
        default=None,
        help=(
            "Use only the first N structures of each length. "
            "Useful for testing or constructing a cluster-growth curve."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output directory for the same run.",
    )

    # Internal worker arguments.
    # These are used by the parent process and should not normally be supplied
    # manually.
    parser.add_argument(
        "--_worker",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--_worker-length",
        type=int,
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
    """
    Find all .cif.gz structures for one peptide length.

    Files are sorted numerically by:
        1. batch number
        2. model number

    This gives a deterministic input order.
    """

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

    designs.sort(key=lambda design: design.sort_key)

    return designs


def check_for_duplicate_ids(
    designs: list[DesignFile],
) -> None:
    """
    Make sure no two files have the same (batch, model) identifier.
    """

    seen: set[tuple[int, int]] = set()

    for design in designs:
        key = (design.batch, design.model)

        if key in seen:
            raise ValueError(
                "Duplicate RFD3 batch/model identifier found: "
                f"batch={design.batch}, model={design.model}"
            )

        seen.add(key)


def write_input_manifest(
    designs: list[DesignFile],
    output_path: Path,
) -> None:
    """
    Save the exact structures and ordering used in the clustering run.
    """

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow(
            [
                "input_rank",
                "peptide_length",
                "batch",
                "model",
                "filename",
                "original_path",
                "external_score",
            ]
        )

        for rank, design in enumerate(designs):
            writer.writerow(
                [
                    rank,
                    design.length,
                    design.batch,
                    design.model,
                    design.path.name,
                    str(design.path),
                    float(rank),
                ]
            )


# =============================================================================
# Temporary Rosetta input preparation
# =============================================================================


def decompress_designs(
    designs: list[DesignFile],
    scratch_dir: Path,
) -> list[Path]:
    """
    Decompress .cif.gz files into temporary .cif files.

    We do this explicitly rather than relying on whether a particular
    PyRosetta build supports gzip-compressed mmCIF input directly.
    """

    output_paths: list[Path] = []
    total = len(designs)

    print(
        f"Decompressing {total:,} mmCIF files...",
        flush=True,
    )

    for index, design in enumerate(designs, start=1):
        output_name = design.path.name.removesuffix(".gz")
        output_path = scratch_dir / output_name

        with gzip.open(design.path, "rb") as source:
            with output_path.open("wb") as destination:
                shutil.copyfileobj(source, destination)

        output_paths.append(output_path)

        if index % 1000 == 0 or index == total:
            print(
                f"  {index:,}/{total:,}",
                flush=True,
            )

    return output_paths


def write_rosetta_input_files(
    cif_paths: list[Path],
    scratch_dir: Path,
) -> tuple[Path, Path]:
    """
    Write the files used by Rosetta.

    rosetta_inputs.txt
        One mmCIF filename per line.

    alternative_scores.txt
        One filename and one score per line.

    The score is simply the deterministic input rank:

        first structure  -> 0
        second structure -> 1
        third structure  -> 2
        ...

    EnergyBasedClusteringProtocol uses score order when selecting new cluster
    centers. These artificial scores therefore make the greedy clustering
    order deterministic without calculating Rosetta energies.
    """

    input_list_path = scratch_dir / "rosetta_inputs.txt"
    score_file_path = scratch_dir / "alternative_scores.txt"

    with input_list_path.open("w") as input_handle:
        with score_file_path.open("w") as score_handle:
            for rank, cif_path in enumerate(cif_paths):
                input_handle.write(f"{cif_path.name}\n")
                score_handle.write(
                    f"{cif_path.name} {float(rank)}\n"
                )

    return input_list_path, score_file_path


# =============================================================================
# PyRosetta clustering
# =============================================================================


def run_rosetta_clustering(
    *,
    input_list_path: Path,
    score_file_path: Path,
    scratch_dir: Path,
    radius: float,
) -> tuple[int, str, float]:
    """
    Run Rosetta's EnergyBasedClusteringProtocol.

    The actual structural comparison, cyclic-permutation testing, and greedy
    clustering are performed by Rosetta.

    Rosetta output is intentionally limited to ONE structure from ONE cluster.
    The limit applies only to disk output; Rosetta still generates the full
    set of clusters internally.
    """

    import pyrosetta
    from pyrosetta.rosetta.protocols import energy_based_clustering as ebc

    original_working_directory = Path.cwd()

    try:
        # The Rosetta input list contains basenames, so work from the directory
        # that contains the temporary mmCIF files.
        os.chdir(scratch_dir)

        pyrosetta.init(
            " ".join(
                [
                    f"-in:file:l {input_list_path.name}",
                    "-in:file:fullatom",
                    "-mute all",
                ]
            )
        )

        # False means that this object will not read clustering settings from
        # Rosetta's global command-line options. We set the settings explicitly.
        options = ebc.EnergyBasedClusteringOptions(False)

        # ---------------------------------------------------------------------
        # Distance metric
        # ---------------------------------------------------------------------

        options.cluster_by_ = ebc.EBC_bb_cartesian
        options.cluster_radius_ = float(radius)

        # Do not add C-beta atoms to the backbone comparison.
        options.use_CB_ = False

        # ---------------------------------------------------------------------
        # Cyclic peptide handling
        # ---------------------------------------------------------------------

        options.cyclic_ = True
        options.cluster_cyclic_permutations_ = True

        # Test every one-residue cyclic shift:
        #
        #   1 2 3 4 5
        #   2 3 4 5 1
        #   3 4 5 1 2
        #   ...
        #
        options.cyclic_permutation_offset_ = 1

        # ---------------------------------------------------------------------
        # Do not modify the generated backbones
        # ---------------------------------------------------------------------

        options.prerelax_ = False
        options.mutate_to_ala_ = False

        # ---------------------------------------------------------------------
        # Deterministic cluster-center ordering
        # ---------------------------------------------------------------------

        options.path_to_scores_file_ = str(
            score_file_path.resolve()
        )

        # ---------------------------------------------------------------------
        # IMPORTANT: minimize Rosetta's PDB-output path
        # ---------------------------------------------------------------------
        #
        # These settings DO NOT limit the actual clustering calculation.
        #
        # Rosetta's documentation defines these as limits on what gets WRITTEN
        # TO DISK. More structures can belong to a cluster and more clusters can
        # be generated internally.
        #
        # We only need n_clusters_from_last_run() for the diversity analysis.
        #
        # Limiting output also avoids Rosetta trying to reconstruct and align
        # every cyclically permuted member, which is the code path that triggered
        # the rms_util.cc numerical assertion in this PyRosetta build.
        #
        options.limit_structures_per_cluster_ = 1
        options.limit_clusters_ = 1

        options.silent_output_ = False

        # Put any diagnostic Rosetta output inside the temporary directory.
        # It will be discarded automatically after the run.
        options.output_prefix_ = str(
            scratch_dir.resolve() / "rosetta_"
        )

        protocol = ebc.EnergyBasedClusteringProtocol(options)

        print(
            "Running Rosetta clustering...",
            flush=True,
        )
        print(
            "Rosetta tracer output is muted; no progress lines are expected.",
            flush=True,
        )

        start_time = time.time()

        protocol.go()

        elapsed_seconds = time.time() - start_time

        n_clusters = int(
            protocol.n_clusters_from_last_run()
        )

        if n_clusters <= 0:
            raise RuntimeError(
                "Rosetta completed but reported zero clusters."
            )

        try:
            pyrosetta_version = str(pyrosetta.version())
        except Exception:
            pyrosetta_version = "unknown"

        print(
            (
                f"Rosetta finished: {n_clusters:,} clusters "
                f"in {elapsed_seconds / 60:.2f} minutes."
            ),
            flush=True,
        )

        return (
            n_clusters,
            pyrosetta_version,
            elapsed_seconds,
        )

    finally:
        os.chdir(original_working_directory)


# =============================================================================
# Output files
# =============================================================================


def write_result_csv(
    *,
    output_path: Path,
    peptide_length: int,
    n_structures: int,
    n_clusters: int,
    radius: float,
    elapsed_seconds: float,
) -> None:
    """
    Write the main numerical result in a simple plotting-friendly CSV.
    """

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "peptide_length",
                "n_structures",
                "n_clusters",
                "unique_cluster_fraction",
                "cluster_radius_angstrom",
                "elapsed_seconds",
            ],
        )

        writer.writeheader()

        writer.writerow(
            {
                "peptide_length": peptide_length,
                "n_structures": n_structures,
                "n_clusters": n_clusters,
                "unique_cluster_fraction": (
                    n_clusters / n_structures
                ),
                "cluster_radius_angstrom": radius,
                "elapsed_seconds": elapsed_seconds,
            }
        )


def write_summary_json(
    *,
    output_path: Path,
    peptide_length: int,
    n_structures: int,
    n_clusters: int,
    radius: float,
    elapsed_seconds: float,
    pyrosetta_version: str,
    input_dir: Path,
) -> None:
    """
    Write complete run provenance and clustering settings.
    """

    summary = {
        "peptide_length": peptide_length,
        "n_structures": n_structures,
        "n_clusters": n_clusters,
        "unique_cluster_fraction": (
            n_clusters / n_structures
        ),
        "cluster_radius_angstrom": radius,
        "elapsed_seconds": elapsed_seconds,
        "input_directory": str(input_dir.resolve()),
        "method": {
            "software": (
                "PyRosetta EnergyBasedClusteringProtocol"
            ),
            "pyrosetta_version": pyrosetta_version,
            "cluster_by": "bb_cartesian",
            "use_CB": False,
            "cyclic": True,
            "cluster_cyclic_permutations": True,
            "cyclic_permutation_offset": 1,
            "prerelax": False,
            "mutate_to_ala": False,
            "seed_order": (
                "deterministic input order via external scores"
            ),
            "input_sort_order": [
                "batch",
                "model",
            ],
            "rosetta_output": {
                "limit_structures_per_cluster": 1,
                "limit_clusters": 1,
                "reason": (
                    "Disk output is deliberately minimized because cyclic "
                    "member reconstruction/alignment triggered a numerical "
                    "assertion in Rosetta rms_util.cc. These limits apply to "
                    "output only; the complete clustering calculation is used "
                    "for n_clusters."
                ),
            },
        },
    }

    with output_path.open("w") as handle:
        json.dump(
            summary,
            handle,
            indent=2,
        )
        handle.write("\n")


# =============================================================================
# Output directory handling
# =============================================================================


def prepare_run_directory(
    run_dir: Path,
    overwrite: bool,
) -> None:
    """Create an empty directory for one clustering run."""

    if run_dir.exists():
        contains_files = any(run_dir.iterdir())

        if contains_files and not overwrite:
            raise FileExistsError(
                "\n".join(
                    [
                        "Output directory already contains files:",
                        f"  {run_dir}",
                        "",
                        "Use --overwrite to replace it.",
                    ]
                )
            )

        if overwrite:
            shutil.rmtree(run_dir)

    run_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


# =============================================================================
# One-length worker process
# =============================================================================


def worker_main(
    args: argparse.Namespace,
) -> None:
    """
    Cluster one peptide length in a fresh PyRosetta process.
    """

    if args._worker_length is None:
        raise ValueError(
            "Worker mode requires --_worker-length."
        )

    peptide_length = args._worker_length
    input_dir = args.input_dir.resolve()
    output_root = args.output_dir.resolve()

    designs = discover_designs(
        input_dir,
        peptide_length,
    )

    if not designs:
        raise FileNotFoundError(
            (
                f"No {peptide_length}-mer .cif.gz files "
                f"were found in:\n{input_dir}"
            )
        )

    check_for_duplicate_ids(designs)

    if args.max_structures is not None:
        if args.max_structures <= 0:
            raise ValueError(
                "--max-structures must be greater than zero."
            )

        designs = designs[: args.max_structures]

    n_structures = len(designs)

    if args.max_structures is None:
        run_name = f"length_{peptide_length}"
    else:
        run_name = (
            f"length_{peptide_length}_n{n_structures}"
        )

    run_dir = output_root / run_name

    prepare_run_directory(
        run_dir,
        args.overwrite,
    )

    print()
    print("=" * 72, flush=True)
    print(
        f"Clustering {peptide_length}-mer macrocycles",
        flush=True,
    )
    print("=" * 72, flush=True)
    print(
        f"Input directory : {input_dir}",
        flush=True,
    )
    print(
        f"Structures      : {n_structures:,}",
        flush=True,
    )
    print(
        f"Cluster radius  : {args.radius:.3f} Å",
        flush=True,
    )
    print(
        f"Output directory: {run_dir}",
        flush=True,
    )
    print(flush=True)

    # Record the exact structures used before starting Rosetta.
    write_input_manifest(
        designs,
        run_dir / "input_manifest.csv",
    )

    scratch_parent = output_root / "_scratch"

    scratch_parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with tempfile.TemporaryDirectory(
        prefix=f"length_{peptide_length}_",
        dir=scratch_parent,
    ) as temporary_directory:

        scratch_dir = Path(
            temporary_directory
        )

        cif_paths = decompress_designs(
            designs,
            scratch_dir,
        )

        (
            input_list_path,
            score_file_path,
        ) = write_rosetta_input_files(
            cif_paths,
            scratch_dir,
        )

        # Save the exact deterministic score ordering permanently.
        shutil.copy2(
            score_file_path,
            run_dir / "alternative_scores.txt",
        )

        (
            n_clusters,
            pyrosetta_version,
            elapsed_seconds,
        ) = run_rosetta_clustering(
            input_list_path=input_list_path,
            score_file_path=score_file_path,
            scratch_dir=scratch_dir,
            radius=args.radius,
        )

    # Temporary mmCIFs and Rosetta PDB output are now deleted.
    # Only the numerical clustering result and provenance are retained.

    write_result_csv(
        output_path=run_dir / "clustering_result.csv",
        peptide_length=peptide_length,
        n_structures=n_structures,
        n_clusters=n_clusters,
        radius=args.radius,
        elapsed_seconds=elapsed_seconds,
    )

    write_summary_json(
        output_path=run_dir / "summary.json",
        peptide_length=peptide_length,
        n_structures=n_structures,
        n_clusters=n_clusters,
        radius=args.radius,
        elapsed_seconds=elapsed_seconds,
        pyrosetta_version=pyrosetta_version,
        input_dir=input_dir,
    )

    try:
        scratch_parent.rmdir()
    except OSError:
        pass

    print()
    print(
        f"Structures      : {n_structures:,}",
        flush=True,
    )
    print(
        f"Clusters        : {n_clusters:,}",
        flush=True,
    )
    print(
        (
            f"Unique fraction : "
            f"{n_clusters / n_structures:.4f}"
        ),
        flush=True,
    )
    print(
        (
            f"Result CSV      : "
            f"{run_dir / 'clustering_result.csv'}"
        ),
        flush=True,
    )
    print(
        (
            f"Summary JSON    : "
            f"{run_dir / 'summary.json'}"
        ),
        flush=True,
    )


# =============================================================================
# Parent process
# =============================================================================


def parent_main(
    args: argparse.Namespace,
) -> None:
    """
    Run each requested peptide length in its own Python process.
    """

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist:\n{input_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for peptide_length in args.lengths:
        command = [
            sys.executable,
            str(SCRIPT_PATH),
            "--_worker",
            "--_worker-length",
            str(peptide_length),
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--radius",
            str(args.radius),
        ]

        if args.max_structures is not None:
            command.extend(
                [
                    "--max-structures",
                    str(args.max_structures),
                ]
            )

        if args.overwrite:
            command.append(
                "--overwrite"
            )

        subprocess.run(
            command,
            check=True,
        )

    print()
    print("=" * 72, flush=True)
    print(
        "All requested peptide lengths complete.",
        flush=True,
    )
    print(
        f"Results: {output_dir}",
        flush=True,
    )
    print("=" * 72, flush=True)


# =============================================================================
# Entry point
# =============================================================================


def main() -> None:
    args = parse_arguments()

    if args.radius <= 0:
        raise ValueError(
            "--radius must be greater than zero."
        )

    if args._worker:
        worker_main(args)
    else:
        parent_main(args)


if __name__ == "__main__":
    main()