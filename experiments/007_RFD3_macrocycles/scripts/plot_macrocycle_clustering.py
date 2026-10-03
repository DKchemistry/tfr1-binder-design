#!/usr/bin/env python3

"""
Plot macrocycle structural diversity as a function of sampling depth.

Reads:

    clustering/length_10/cluster_growth.csv
    clustering/length_12/cluster_growth.csv

and creates an RFpeptides Fig. S2H-style plot of:

    sampled structures
            vs.
    number of structural clusters

Outputs:

    clustering/figures/macrocycle_cluster_growth.pdf
    clustering/figures/macrocycle_cluster_growth.png
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

# Always render without opening a GUI.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter


# =============================================================================
# Paths
# =============================================================================

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1]

DEFAULT_CLUSTERING_DIR = (
    EXPERIMENT_DIR
    / "test_2"
    / "macrocycle_monomer_10k"
    / "clustering"
)

DEFAULT_OUTPUT_DIR = (
    DEFAULT_CLUSTERING_DIR
    / "figures"
)


# =============================================================================
# Data structure
# =============================================================================


@dataclass(frozen=True)
class GrowthPoint:
    peptide_length: int
    n_structures: int
    n_clusters: int


# =============================================================================
# Arguments
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot macrocycle backbone cluster growth."
    )

    parser.add_argument(
        "--clustering-dir",
        type=Path,
        default=DEFAULT_CLUSTERING_DIR,
        help=f"Clustering results directory. Default: {DEFAULT_CLUSTERING_DIR}",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Figure output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=[10, 12],
        help="Peptide lengths to plot. Default: 10 12",
    )

    return parser.parse_args()


# =============================================================================
# Data loading
# =============================================================================


def read_growth_curve(
    clustering_dir: Path,
    peptide_length: int,
) -> list[GrowthPoint]:
    """Read cluster-growth data for one peptide length."""

    input_path = (
        clustering_dir
        / f"length_{peptide_length}"
        / "cluster_growth.csv"
    )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Could not find cluster-growth data:\n{input_path}"
        )

    points: list[GrowthPoint] = []

    with input_path.open() as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            points.append(
                GrowthPoint(
                    peptide_length=int(row["peptide_length"]),
                    n_structures=int(row["n_structures"]),
                    n_clusters=int(row["n_clusters"]),
                )
            )

    if not points:
        raise ValueError(
            f"No data rows found in {input_path}"
        )

    points.sort(key=lambda point: point.n_structures)

    return points


# =============================================================================
# Plotting
# =============================================================================


def make_plot(
    curves: dict[int, list[GrowthPoint]],
    output_dir: Path,
) -> None:
    """Create a simple academic-style cluster-growth plot."""

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(
        figsize=(3.3, 3.1)
    )

    # Different markers make the curves distinguishable even if the figure is
    # viewed or printed in grayscale.
    markers = ["o", "s", "^", "D"]

    for index, peptide_length in enumerate(sorted(curves)):
        points = curves[peptide_length]

        x = [
            point.n_structures
            for point in points
        ]

        y = [
            point.n_clusters
            for point in points
        ]

        ax.plot(
            x,
            y,
            marker=markers[index % len(markers)],
            linewidth=1.5,
            markersize=5,
            label=f"{peptide_length} residues",
        )

    ax.set_xlabel("Sampled structures")
    ax.set_ylabel("Clusters")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        direction="out",
        length=3,
    )

    ax.legend(
        frameon=False,
        loc="best",
    )

    # Show the actual sampling depths used by the analysis.
    sample_sizes = sorted(
        {
            point.n_structures
            for points in curves.values()
            for point in points
        }
    )

    ax.set_xticks(sample_sizes)

    ax.set_xticklabels(
        [f"{value:,}" for value in sample_sizes],
        rotation=30,
        ha="right",
    )

    ax.yaxis.set_major_formatter(
        StrMethodFormatter("{x:,.0f}")
    )

    # As in the RFpeptides panel, the sampling axis starts at zero even though
    # the first measured point is later.
    ax.set_xlim(
        left=0,
        right=max(sample_sizes) * 1.04,
    )

    fig.tight_layout()

    pdf_path = (
        output_dir
        / "macrocycle_cluster_growth.pdf"
    )

    png_path = (
        output_dir
        / "macrocycle_cluster_growth.png"
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    args = parse_arguments()

    clustering_dir = args.clustering_dir.resolve()

    curves = {
        peptide_length: read_growth_curve(
            clustering_dir,
            peptide_length,
        )
        for peptide_length in args.lengths
    }

    print("Cluster-growth data:")

    for peptide_length in sorted(curves):
        print(f"  {peptide_length}-mer")

        for point in curves[peptide_length]:
            print(
                f"    {point.n_structures:>10,} structures -> "
                f"{point.n_clusters:>6,} clusters"
            )

    make_plot(
        curves,
        args.output_dir.resolve(),
    )


if __name__ == "__main__":
    main()