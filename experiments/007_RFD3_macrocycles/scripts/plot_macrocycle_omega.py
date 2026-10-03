#!/usr/bin/env python3

"""
Plot omega-angle histograms for generated macrocycle peptide bonds.

This is an all-design geometry diagnostic. It is intentionally labelled in
peptide-bond observations, not clusters: reproducing RFpeptides Figure S2J
exactly requires one representative from every backbone cluster. Pass the
saved cluster-center manifests to select those structures from the existing
all-design torsion table.

Input:

    ramachandran/all_structures_torsions.csv

Outputs:

    ramachandran/figures/macrocycle_omega.pdf
    ramachandran/figures/macrocycle_omega.png
    ramachandran/figures/macrocycle_cluster_center_omega.pdf
    ramachandran/figures/macrocycle_cluster_center_omega.png
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# Paths and colors
# =============================================================================

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1]

DEFAULT_RAMACHANDRAN_DIR = (
    EXPERIMENT_DIR
    / "test_2"
    / "macrocycle_monomer_10k"
    / "ramachandran"
)

DEFAULT_INPUT_PATH = DEFAULT_RAMACHANDRAN_DIR / "all_structures_torsions.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_RAMACHANDRAN_DIR / "figures"

LENGTH_COLORS = {
    10: "#3f007d",
    12: "#149c7e",
}


# =============================================================================
# Arguments and data loading
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot macrocycle omega-angle histograms."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help=f"Torsion CSV. Default: {DEFAULT_INPUT_PATH}",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for figures. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=[10, 12],
        help="Peptide lengths to plot. Default: 10 12",
    )

    parser.add_argument(
        "--bin-width",
        type=float,
        default=2.0,
        help="Histogram bin width in degrees. Default: 2",
    )

    parser.add_argument(
        "--cluster-center-manifests",
        type=Path,
        nargs="+",
        default=None,
        help=(
            "Optional cluster-center manifest.csv files. When supplied, "
            "only torsions from those representative structures are plotted."
        ),
    )

    return parser.parse_args()


def read_omega_angles(
    input_path: Path,
    requested_lengths: set[int],
    selected_filenames: dict[int, set[str]] | None,
) -> dict[int, list[float]]:
    """Read omega angles, grouped by peptide length."""

    angles = {length: [] for length in requested_lengths}

    with input_path.open() as handle:
        reader = csv.DictReader(handle)

        required_columns = {"peptide_length", "filename", "omega_deg"}
        missing_columns = required_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            peptide_length = int(row["peptide_length"])

            if peptide_length not in angles:
                continue

            if (
                selected_filenames is not None
                and row["filename"]
                not in selected_filenames.get(peptide_length, set())
            ):
                continue

            angles[peptide_length].append(float(row["omega_deg"]))

    for peptide_length, values in angles.items():
        if not values:
            raise ValueError(
                f"No omega angles found for {peptide_length}-mers."
            )

    return angles


def read_cluster_center_manifests(
    manifest_paths: list[Path],
    requested_lengths: set[int],
) -> dict[int, set[str]]:
    """Read the representative filenames selected by clustering."""

    selected_filenames = {
        peptide_length: set()
        for peptide_length in requested_lengths
    }

    for manifest_path in manifest_paths:
        with manifest_path.open() as handle:
            reader = csv.DictReader(handle)
            required_columns = {
                "cluster_id",
                "peptide_length",
                "filename",
            }
            missing_columns = required_columns - set(reader.fieldnames or [])

            if missing_columns:
                raise ValueError(
                    f"{manifest_path} is missing required columns: "
                    + ", ".join(sorted(missing_columns))
                )

            for row in reader:
                peptide_length = int(row["peptide_length"])

                if peptide_length not in requested_lengths:
                    continue

                filename = row["filename"]

                if filename in selected_filenames[peptide_length]:
                    raise ValueError(
                        f"Duplicate {peptide_length}-mer representative "
                        f"{filename} in the supplied manifests."
                    )

                selected_filenames[peptide_length].add(filename)

    for peptide_length, filenames in selected_filenames.items():
        if not filenames:
            raise ValueError(
                f"No {peptide_length}-mer representatives were found in "
                "the supplied manifests."
            )

    return selected_filenames


# =============================================================================
# Plotting
# =============================================================================


def configure_matplotlib() -> None:
    """Apply the same compact style as the Ramachandran figures."""

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_omega(
    angles_by_length: dict[int, list[float]],
    peptide_lengths: list[int],
    bin_width: float,
    cluster_counts: dict[int, int] | None,
) -> plt.Figure:
    """Create linear- and logarithmic-scale omega histograms."""

    bin_edges = np.arange(
        -180.0,
        180.0 + bin_width,
        bin_width,
    )

    figure, axes = plt.subplots(
        nrows=len(peptide_lengths),
        ncols=2,
        figsize=(8.0, 2.8 * len(peptide_lengths)),
        sharex=True,
        squeeze=False,
    )

    for row_index, peptide_length in enumerate(peptide_lengths):
        angles = np.asarray(angles_by_length[peptide_length])
        weights = np.full(len(angles), 100.0 / len(angles))
        color = LENGTH_COLORS.get(peptide_length, "#333333")

        for column_index, use_log_scale in enumerate((False, True)):
            axis = axes[row_index, column_index]

            axis.hist(
                angles,
                bins=bin_edges,
                weights=weights,
                color=color,
                edgecolor="none",
            )

            axis.axvspan(-30, 30, color="#f2b134", alpha=0.12)
            axis.axvline(-150, color="#777777", linewidth=0.7, linestyle="--")
            axis.axvline(150, color="#777777", linewidth=0.7, linestyle="--")
            axis.set_xlim(-180, 180)
            axis.set_xticks([-180, -90, 0, 90, 180])
            axis.tick_params(direction="out", length=3)

            if use_log_scale:
                axis.set_yscale("log")
                axis.set_title("Log scale")
            else:
                axis.set_title("Linear scale")

            if column_index == 0:
                observation_label = (
                    "Representative bonds"
                    if cluster_counts is not None
                    else "Peptide bonds"
                )
                axis.set_ylabel(
                    f"{peptide_length}-mer\n"
                    f"{observation_label} per bin (%)"
                )

            if row_index == len(peptide_lengths) - 1:
                axis.set_xlabel(r"$\omega$ (degrees)")

            annotation = f"n={len(angles):,} bonds"

            if cluster_counts is not None:
                annotation = (
                    f"{cluster_counts[peptide_length]:,} clusters\n"
                    f"{len(angles):,} bonds"
                )

            axis.text(
                0.03,
                0.92,
                annotation,
                transform=axis.transAxes,
                va="top",
                ha="left",
                fontsize=8,
            )

    if cluster_counts is not None:
        figure.suptitle(
            "Omega angles from backbone-cluster representatives",
            fontsize=11,
        )
        figure.tight_layout(rect=(0, 0, 1, 0.97))
    else:
        figure.tight_layout()

    return figure


def save_figure(
    figure: plt.Figure,
    output_dir: Path,
    stem: str,
) -> None:
    """Save the figure as PDF and high-resolution PNG."""

    output_dir.mkdir(parents=True, exist_ok=True)

    for suffix in ("pdf", "png"):
        output_path = output_dir / f"{stem}.{suffix}"
        save_options = {"bbox_inches": "tight"}

        if suffix == "png":
            save_options["dpi"] = 300

        figure.savefig(output_path, **save_options)
        print(f"Saved: {output_path}")

    plt.close(figure)


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    args = parse_arguments()

    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()
    peptide_lengths = sorted(set(args.lengths))

    if not input_path.is_file():
        raise FileNotFoundError(f"Torsion CSV does not exist:\n{input_path}")

    if args.bin_width <= 0 or 360.0 % args.bin_width != 0:
        raise ValueError(
            "--bin-width must be positive and divide 360 degrees evenly."
        )

    selected_filenames = None

    if args.cluster_center_manifests is not None:
        manifest_paths = [
            path.resolve()
            for path in args.cluster_center_manifests
        ]

        for manifest_path in manifest_paths:
            if not manifest_path.is_file():
                raise FileNotFoundError(
                    f"Cluster-center manifest does not exist:\n{manifest_path}"
                )

        selected_filenames = read_cluster_center_manifests(
            manifest_paths,
            set(peptide_lengths),
        )

    angles_by_length = read_omega_angles(
        input_path,
        set(peptide_lengths),
        selected_filenames,
    )

    cluster_counts = None
    output_stem = "macrocycle_omega"

    if selected_filenames is not None:
        cluster_counts = {
            peptide_length: len(filenames)
            for peptide_length, filenames in selected_filenames.items()
        }
        output_stem = "macrocycle_cluster_center_omega"

        for peptide_length in peptide_lengths:
            expected_angle_count = (
                cluster_counts[peptide_length]
                * peptide_length
            )
            observed_angle_count = len(angles_by_length[peptide_length])

            if observed_angle_count != expected_angle_count:
                raise ValueError(
                    f"Expected {expected_angle_count:,} omega angles from "
                    f"the {cluster_counts[peptide_length]:,} saved "
                    f"{peptide_length}-mer centers, but found "
                    f"{observed_angle_count:,}."
                )

    configure_matplotlib()
    figure = plot_omega(
        angles_by_length,
        peptide_lengths,
        args.bin_width,
        cluster_counts,
    )
    save_figure(figure, output_dir, output_stem)


if __name__ == "__main__":
    main()
