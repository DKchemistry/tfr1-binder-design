#!/usr/bin/env python3

"""Plot the RFD3 -> LigandMPNN -> AfCyc monomer analysis.

Reads the tab-separated tables written by
``analyze_afcyc_monomer_self_consistency.py`` and creates three figures:

* AfCyc pLDDT versus backbone self-consistency RMSD;
* backbone success rate as a function of peptide length;
* backbone success rate as a function of LigandMPNN sequence attempts.

The success definition follows RFpeptides Figure 1c: a backbone succeeds when
at least one sequence has pLDDT > 0.8 and backbone RMSD < 2.0 Angstrom.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

# Always render without opening a graphical window.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


# =============================================================================
# Paths
# =============================================================================

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1] / "monomer_self_consistency"
DEFAULT_ANALYSIS_DIR = EXPERIMENT_DIR / "analysis" / "run_01"
DEFAULT_OUTPUT_DIR = DEFAULT_ANALYSIS_DIR / "figures"


# =============================================================================
# Arguments and data loading
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=DEFAULT_ANALYSIS_DIR,
        help=f"Analysis table directory. Default: {DEFAULT_ANALYSIS_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Figure output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated table and return its rows."""

    if not path.is_file():
        raise FileNotFoundError(f"Could not find analysis table:\n{path}")

    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        raise ValueError(f"No data rows found in {path}")

    return rows


def read_cutoffs(summary_path: Path) -> tuple[float, float]:
    """Read the pLDDT and RMSD cutoffs used by the analysis."""

    if not summary_path.is_file():
        raise FileNotFoundError(f"Could not find analysis summary:\n{summary_path}")

    with summary_path.open() as handle:
        summary = json.load(handle)

    method = summary["method"]
    return float(method["plddt_cutoff"]), float(method["rmsd_cutoff_angstrom"])


# =============================================================================
# Shared plotting style
# =============================================================================


def set_plot_style() -> None:
    """Use the same compact style as the other experiment 007 plots."""

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


def finish_axes(ax: plt.Axes) -> None:
    """Apply shared axis formatting."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3)


def save_figure(fig: plt.Figure, output_dir: Path, filename_stem: str) -> None:
    """Save one figure as both PDF and high-resolution PNG."""

    output_dir.mkdir(parents=True, exist_ok=True)

    pdf_path = output_dir / f"{filename_stem}.pdf"
    png_path = output_dir / f"{filename_stem}.png"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


# =============================================================================
# Figures
# =============================================================================


def plot_plddt_vs_rmsd(
    prediction_rows: list[dict[str, str]],
    plddt_cutoff: float,
    rmsd_cutoff: float,
    output_dir: Path,
) -> None:
    """Plot every completed AfCyc sequence prediction."""

    rows_by_length: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in prediction_rows:
        rows_by_length[int(row["peptide_length"])].append(row)

    fig, ax = plt.subplots(figsize=(4.1, 3.3))

    for peptide_length in sorted(rows_by_length):
        rows = rows_by_length[peptide_length]
        rmsd_values = [float(row["backbone_rmsd_angstrom"]) for row in rows]
        plddt_values = [float(row["plddt"]) for row in rows]

        ax.scatter(
            rmsd_values,
            plddt_values,
            s=9,
            alpha=0.35,
            edgecolors="none",
            label=f"{peptide_length}-mer",
        )

    # Dashed lines show the two parts of the RFpeptides success definition.
    ax.axvline(rmsd_cutoff, color="0.35", linestyle="--", linewidth=1.0)
    ax.axhline(plddt_cutoff, color="0.35", linestyle="--", linewidth=1.0)

    ax.set_xlabel(r"Backbone scRMSD ($\AA$)")
    ax.set_ylabel("AfCyc pLDDT")
    ax.set_xlim(left=0)
    ax.set_ylim(0.55, 1.0)
    finish_axes(ax)

    ax.legend(
        frameon=False,
        loc="lower right",
        ncol=2,
        markerscale=1.8,
    )

    fig.tight_layout()
    save_figure(fig, output_dir, "afcyc_plddt_vs_backbone_rmsd")


def plot_success_by_length(
    length_rows: list[dict[str, str]],
    output_dir: Path,
) -> None:
    """Plot the RFpeptides-style backbone success rate by peptide length."""

    complete_rows = [row for row in length_rows if row["campaign_complete"] == "True"]
    if len(complete_rows) != len(length_rows):
        raise ValueError(
            "success_by_length.tsv contains an incomplete campaign; rerun the "
            "analysis after AfCyc finishes before making the final plot"
        )

    complete_rows.sort(key=lambda row: int(row["peptide_length"]))
    peptide_lengths = [int(row["peptide_length"]) for row in complete_rows]
    success_rates = [
        float(row["rfpeptides_style_success_fraction"])
        for row in complete_rows
    ]

    fig, ax = plt.subplots(figsize=(3.3, 3.1))

    ax.plot(
        peptide_lengths,
        success_rates,
        marker="o",
        linewidth=1.5,
        markersize=5,
    )

    ax.set_xlabel("Peptide length")
    ax.set_ylabel("Successful backbones")
    ax.set_xticks(peptide_lengths)
    ax.set_ylim(0.0, 1.02)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    finish_axes(ax)

    fig.tight_layout()
    save_figure(fig, output_dir, "afcyc_success_by_length")


def plot_success_by_attempts(
    attempt_rows: list[dict[str, str]],
    output_dir: Path,
) -> None:
    """Plot success after the first one through eight sequence attempts."""

    rows_by_length: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in attempt_rows:
        rows_by_length[int(row["peptide_length"])].append(row)

    fig, ax = plt.subplots(figsize=(4.1, 3.3))
    markers = ["o", "s", "^", "D", "v", "P"]

    for index, peptide_length in enumerate(sorted(rows_by_length)):
        rows = rows_by_length[peptide_length]
        rows.sort(key=lambda row: int(row["sequence_attempts"]))

        attempts = [int(row["sequence_attempts"]) for row in rows]
        success_rates = [float(row["success_fraction"]) for row in rows]

        ax.plot(
            attempts,
            success_rates,
            marker=markers[index % len(markers)],
            linewidth=1.5,
            markersize=4,
            label=f"{peptide_length}-mer",
        )

    ax.set_xlabel("LigandMPNN sequence attempts")
    ax.set_ylabel("Successful backbones")
    ax.set_xticks(range(1, 9))
    ax.set_ylim(0.0, 1.02)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    finish_axes(ax)

    ax.legend(
        frameon=False,
        loc="lower right",
        ncol=2,
    )

    fig.tight_layout()
    save_figure(fig, output_dir, "afcyc_success_by_sequence_attempts")


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    args = parse_arguments()
    analysis_dir = args.analysis_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    prediction_rows = read_tsv(analysis_dir / "prediction_results.tsv")
    length_rows = read_tsv(analysis_dir / "success_by_length.tsv")
    attempt_rows = read_tsv(analysis_dir / "success_by_sequence_attempts.tsv")
    plddt_cutoff, rmsd_cutoff = read_cutoffs(
        analysis_dir / "analysis_summary.json"
    )

    set_plot_style()
    plot_plddt_vs_rmsd(
        prediction_rows,
        plddt_cutoff,
        rmsd_cutoff,
        output_dir,
    )
    plot_success_by_length(length_rows, output_dir)
    plot_success_by_attempts(attempt_rows, output_dir)


if __name__ == "__main__":
    main()
