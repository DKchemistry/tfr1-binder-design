#!/usr/bin/env python3

"""
Plot Ramachandran distributions for RFD3-generated macrocycles.

The primary figure mirrors RFpeptides Supplementary Figure S2I:

    - one panel for 10-residue macrocycles
    - one panel for 12-residue macrocycles
    - every residue from every generated structure

A second figure separates the same observations into the six residue classes
used by MolProbity's Top8000 reference distributions.

Inputs:

    ramachandran/all_structures_torsions.csv

Outputs:

    ramachandran/figures/macrocycle_ramachandran.pdf
    ramachandran/figures/macrocycle_ramachandran.png
    ramachandran/figures/macrocycle_ramachandran_by_residue_class.pdf
    ramachandran/figures/macrocycle_ramachandran_by_residue_class.png
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
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np


# =============================================================================
# Paths and plotting constants
# =============================================================================

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1]

DEFAULT_RAMACHANDRAN_DIR = (
    EXPERIMENT_DIR
    / "test_2"
    / "macrocycle_monomer_10k"
    / "ramachandran"
)

DEFAULT_INPUT_PATH = (
    DEFAULT_RAMACHANDRAN_DIR
    / "all_structures_torsions.csv"
)

DEFAULT_OUTPUT_DIR = (
    DEFAULT_RAMACHANDRAN_DIR
    / "figures"
)

ANGLE_LIMITS = (-180.0, 180.0)

RESIDUE_CLASS_ORDER = [
    "general",
    "glycine",
    "cis-proline",
    "trans-proline",
    "pre-proline",
    "isoleucine or valine",
]

RESIDUE_CLASS_LABELS = {
    "general": "General",
    "glycine": "Glycine",
    "cis-proline": "cis-Pro",
    "trans-proline": "trans-Pro",
    "pre-proline": "Pre-proline",
    "isoleucine or valine": "Ile/Val",
}

LENGTH_COLOR_INDICES = {
    10: 0,
    12: 1,
}


# =============================================================================
# Data structure
# =============================================================================


@dataclass(frozen=True)
class AngleObservation:
    """One residue's phi/psi observation."""

    peptide_length: int
    residue_class: str
    phi_deg: float
    psi_deg: float


# =============================================================================
# Arguments
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot cyclic macrocycle Ramachandran distributions."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help=f"Torsion CSV from the calculation script. Default: {DEFAULT_INPUT_PATH}",
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
        "--bins",
        type=int,
        default=120,
        help="Number of equal-width bins per angle axis. Default: 120",
    )

    return parser.parse_args()


# =============================================================================
# Data loading and binning
# =============================================================================


def read_observations(
    input_path: Path,
    requested_lengths: set[int],
) -> list[AngleObservation]:
    """Read the phi/psi columns needed for plotting."""

    observations: list[AngleObservation] = []

    with input_path.open() as handle:
        reader = csv.DictReader(handle)

        required_columns = {
            "peptide_length",
            "ramachandran_class",
            "phi_deg",
            "psi_deg",
        }

        missing_columns = required_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            peptide_length = int(row["peptide_length"])

            if peptide_length not in requested_lengths:
                continue

            observations.append(
                AngleObservation(
                    peptide_length=peptide_length,
                    residue_class=row["ramachandran_class"],
                    phi_deg=float(row["phi_deg"]),
                    psi_deg=float(row["psi_deg"]),
                )
            )

    return observations


def select_angles(
    observations: list[AngleObservation],
    *,
    peptide_length: int,
    residue_class: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Select phi and psi arrays for one requested population."""

    selected = [
        observation
        for observation in observations
        if observation.peptide_length == peptide_length
        and (
            residue_class is None
            or observation.residue_class == residue_class
        )
    ]

    if not selected:
        population = f"{peptide_length}-mer"

        if residue_class is not None:
            population += f" {residue_class}"

        raise ValueError(f"No observations found for {population}.")

    phi = np.asarray(
        [observation.phi_deg for observation in selected],
        dtype=float,
    )

    psi = np.asarray(
        [observation.psi_deg for observation in selected],
        dtype=float,
    )

    return phi, psi


def probability_density(
    phi: np.ndarray,
    psi: np.ndarray,
    bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate a normalized two-dimensional histogram."""

    density, phi_edges, psi_edges = np.histogram2d(
        phi,
        psi,
        bins=bins,
        range=[ANGLE_LIMITS, ANGLE_LIMITS],
        density=True,
    )

    return density, phi_edges, psi_edges


def make_length_colormap(peptide_length: int) -> LinearSegmentedColormap:
    """Create a white-to-color map for one peptide length."""

    default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    color_index = LENGTH_COLOR_INDICES.get(peptide_length, 0)
    final_color = default_colors[color_index % len(default_colors)]

    return LinearSegmentedColormap.from_list(
        f"length_{peptide_length}",
        ["#ffffff", final_color],
    )


def common_log_norm(
    densities: list[np.ndarray],
) -> LogNorm:
    """Use one logarithmic density scale across related panels."""

    positive_values = np.concatenate(
        [density[density > 0] for density in densities]
    )

    maximum = float(positive_values.max())
    minimum = max(
        float(positive_values.min()),
        maximum * 1e-4,
    )

    return LogNorm(
        vmin=minimum,
        vmax=maximum,
        clip=True,
    )


# =============================================================================
# Plot formatting
# =============================================================================


def configure_matplotlib() -> None:
    """Apply consistent academic-style formatting."""

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def format_ramachandran_axis(
    axis: plt.Axes,
    *,
    show_x_label: bool,
    show_y_label: bool,
) -> None:
    """Format one Ramachandran axis."""

    ticks = [-180, -90, 0, 90, 180]

    axis.set_xlim(ANGLE_LIMITS)
    axis.set_ylim(ANGLE_LIMITS)
    axis.set_xticks(ticks)
    axis.set_yticks(ticks)
    axis.set_aspect("equal")

    axis.set_xlabel(r"$\Phi$ (degrees)" if show_x_label else "")
    axis.set_ylabel(r"$\Psi$ (degrees)" if show_y_label else "")

    axis.tick_params(
        direction="out",
        length=3,
    )


def draw_density(
    axis: plt.Axes,
    density: np.ndarray,
    phi_edges: np.ndarray,
    psi_edges: np.ndarray,
    *,
    peptide_length: int,
    norm: LogNorm,
) -> None:
    """Draw one probability-density panel."""

    masked_density = np.ma.masked_where(
        density <= 0,
        density,
    )

    axis.pcolormesh(
        phi_edges,
        psi_edges,
        masked_density.T,
        cmap=make_length_colormap(peptide_length),
        norm=norm,
        shading="auto",
        rasterized=True,
    )


def save_figure(
    figure: plt.Figure,
    output_dir: Path,
    stem: str,
) -> None:
    """Save one figure as both PDF and high-resolution PNG."""

    pdf_path = output_dir / f"{stem}.pdf"
    png_path = output_dir / f"{stem}.png"

    figure.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    figure.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(figure)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# =============================================================================
# Figures
# =============================================================================


def plot_all_residues(
    observations: list[AngleObservation],
    peptide_lengths: list[int],
    bins: int,
    output_dir: Path,
) -> None:
    """Create the RFpeptides Figure S2I-style plot."""

    populations = {
        peptide_length: select_angles(
            observations,
            peptide_length=peptide_length,
        )
        for peptide_length in peptide_lengths
    }

    histograms = {
        peptide_length: probability_density(phi, psi, bins)
        for peptide_length, (phi, psi) in populations.items()
    }

    norm = common_log_norm(
        [histogram[0] for histogram in histograms.values()]
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=len(peptide_lengths),
        figsize=(4.0 * len(peptide_lengths), 3.6),
        squeeze=False,
    )

    for column_index, peptide_length in enumerate(peptide_lengths):
        axis = axes[0, column_index]
        phi, _ = populations[peptide_length]
        density, phi_edges, psi_edges = histograms[peptide_length]

        draw_density(
            axis,
            density,
            phi_edges,
            psi_edges,
            peptide_length=peptide_length,
            norm=norm,
        )

        axis.set_title(
            f"{peptide_length} residues "
            f"({len(phi):,} residue observations)"
        )

        format_ramachandran_axis(
            axis,
            show_x_label=True,
            show_y_label=(column_index == 0),
        )

    figure.tight_layout()

    save_figure(
        figure,
        output_dir,
        "macrocycle_ramachandran",
    )


def plot_residue_classes(
    observations: list[AngleObservation],
    peptide_lengths: list[int],
    bins: int,
    output_dir: Path,
) -> None:
    """Create a residue-class-aware companion figure."""

    populations: dict[tuple[int, str], tuple[np.ndarray, np.ndarray]] = {}
    histograms: dict[
        tuple[int, str],
        tuple[np.ndarray, np.ndarray, np.ndarray],
    ] = {}

    for peptide_length in peptide_lengths:
        for residue_class in RESIDUE_CLASS_ORDER:
            key = (peptide_length, residue_class)

            populations[key] = select_angles(
                observations,
                peptide_length=peptide_length,
                residue_class=residue_class,
            )

            histograms[key] = probability_density(
                *populations[key],
                bins,
            )

    norm = common_log_norm(
        [histogram[0] for histogram in histograms.values()]
    )

    figure, axes = plt.subplots(
        nrows=len(peptide_lengths),
        ncols=len(RESIDUE_CLASS_ORDER),
        figsize=(17.0, 3.1 * len(peptide_lengths)),
        squeeze=False,
    )

    for row_index, peptide_length in enumerate(peptide_lengths):
        for column_index, residue_class in enumerate(RESIDUE_CLASS_ORDER):
            axis = axes[row_index, column_index]
            key = (peptide_length, residue_class)
            phi, _ = populations[key]
            density, phi_edges, psi_edges = histograms[key]

            draw_density(
                axis,
                density,
                phi_edges,
                psi_edges,
                peptide_length=peptide_length,
                norm=norm,
            )

            if row_index == 0:
                axis.set_title(RESIDUE_CLASS_LABELS[residue_class])

            axis.text(
                0.03,
                0.96,
                f"{peptide_length}-mer\nn={len(phi):,}",
                transform=axis.transAxes,
                va="top",
                ha="left",
                fontsize=8,
            )

            format_ramachandran_axis(
                axis,
                show_x_label=(row_index == len(peptide_lengths) - 1),
                show_y_label=(column_index == 0),
            )

    figure.tight_layout()

    save_figure(
        figure,
        output_dir,
        "macrocycle_ramachandran_by_residue_class",
    )


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    args = parse_arguments()

    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()
    peptide_lengths = sorted(set(args.lengths))

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Torsion CSV does not exist:\n{input_path}"
        )

    if args.bins <= 0:
        raise ValueError("--bins must be greater than zero.")

    observations = read_observations(
        input_path,
        set(peptide_lengths),
    )

    found_lengths = {
        observation.peptide_length
        for observation in observations
    }

    missing_lengths = set(peptide_lengths) - found_lengths

    if missing_lengths:
        raise ValueError(
            "Input CSV has no observations for length(s): "
            + ", ".join(str(length) for length in sorted(missing_lengths))
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    print(f"Input: {input_path}")
    print(f"Observations: {len(observations):,}")
    print(f"Bins per axis: {args.bins}")

    plot_all_residues(
        observations,
        peptide_lengths,
        args.bins,
        output_dir,
    )

    plot_residue_classes(
        observations,
        peptide_lengths,
        args.bins,
        output_dir,
    )


if __name__ == "__main__":
    main()
