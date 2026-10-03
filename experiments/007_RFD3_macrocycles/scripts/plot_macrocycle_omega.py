#!/usr/bin/env python3

"""
Plot omega-angle histograms for all generated macrocycle peptide bonds.

This is an all-design geometry diagnostic. It is intentionally labelled in
peptide-bond observations, not clusters: reproducing RFpeptides Figure S2J
exactly would first require saving one representative from every backbone
cluster.

Input:

    ramachandran/all_structures_torsions.csv

Outputs:

    ramachandran/figures/macrocycle_omega.pdf
    ramachandran/figures/macrocycle_omega.png
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

    return parser.parse_args()


def read_omega_angles(
    input_path: Path,
    requested_lengths: set[int],
) -> dict[int, list[float]]:
    """Read omega angles, grouped by peptide length."""

    angles = {length: [] for length in requested_lengths}

    with input_path.open() as handle:
        reader = csv.DictReader(handle)

        required_columns = {"peptide_length", "omega_deg"}
        missing_columns = required_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            peptide_length = int(row["peptide_length"])

            if peptide_length in angles:
                angles[peptide_length].append(float(row["omega_deg"]))

    for peptide_length, values in angles.items():
        if not values:
            raise ValueError(
                f"No omega angles found for {peptide_length}-mers."
            )

    return angles


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
                axis.set_ylabel(
                    f"{peptide_length}-mer\nPeptide bonds per bin (%)"
                )

            if row_index == len(peptide_lengths) - 1:
                axis.set_xlabel(r"$\omega$ (degrees)")

            axis.text(
                0.03,
                0.92,
                f"n={len(angles):,}",
                transform=axis.transAxes,
                va="top",
                ha="left",
                fontsize=8,
            )

    figure.tight_layout()
    return figure


def save_figure(
    figure: plt.Figure,
    output_dir: Path,
) -> None:
    """Save the figure as PDF and high-resolution PNG."""

    output_dir.mkdir(parents=True, exist_ok=True)

    for suffix in ("pdf", "png"):
        output_path = output_dir / f"macrocycle_omega.{suffix}"
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

    angles_by_length = read_omega_angles(
        input_path,
        set(peptide_lengths),
    )

    configure_matplotlib()
    figure = plot_omega(
        angles_by_length,
        peptide_lengths,
        args.bin_width,
    )
    save_figure(figure, output_dir)


if __name__ == "__main__":
    main()
