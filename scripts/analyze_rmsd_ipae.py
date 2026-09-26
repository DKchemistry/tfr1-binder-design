#!/usr/bin/env python3

"""Analyze target-aligned binder C-alpha RMSD and oracle iPAE values.

For every design found below ``--oracle``, the relaxed target is superposed on
the oracle target with a least-squares Kabsch fit.  RMSD is then calculated
between the relaxed and oracle binder C-alpha atoms.  Residues are paired in
PDB order, which supports the different residue numbering used by the two
stages of this pipeline.  Biotite handles PDB parsing, superposition,
coordinate transformation, and RMSD calculation.

Outputs:
    results.csv             Per-design values, input paths, and errors
    summary.json            Summary statistics and analysis configuration
    rmsd_histogram.png      Binder C-alpha RMSD distribution
    ipae_histogram.png      i_pae_normalized distribution
    rmsd_vs_ipae.png        Relationship between the two metrics
    rmsd_vs_ipae.html       Interactive relationship between the metrics
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go


HOTSPOT_CATEGORIES = ("≤ 4 Å", "> 4 Å")
HOTSPOT_COLORS = {
    "≤ 4 Å": "#1b9e77",
    "> 4 Å": "#d95f02",
}


def parse_hotspot_residues(value: str) -> tuple[int, ...]:
    match = re.fullmatch(r"(-?\d+)(?:-(-?\d+))?", value.strip())
    if match is None:
        raise argparse.ArgumentTypeError(
            "hotspot residues must be a single residue or an inclusive range, "
            "such as 23 or 20-25"
        )

    start = int(match.group(1))
    end = int(match.group(2)) if match.group(2) is not None else start
    if end < start:
        raise argparse.ArgumentTypeError(
            "hotspot residue range must end at or after its start"
        )
    return tuple(range(start, end + 1))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate target-aligned binder C-alpha RMSD for paired relaxed/"
            "oracle structures and plot RMSD and normalized iPAE distributions."
        )
    )
    parser.add_argument(
        "--oracle",
        required=True,
        type=Path,
        help="Oracle output directory containing rfdiffusion_x directories.",
    )
    parser.add_argument(
        "--relaxed",
        required=True,
        type=Path,
        help=(
            "MPNN/relax output root containing rfdiffusion_x directories, or "
            "one relaxed.pdb for a single-design analysis."
        ),
    )
    parser.add_argument(
        "--output-dir",
        "--output_dir",
        dest="output_dir",
        required=True,
        type=Path,
        help="Directory for figures and statistics.",
    )
    parser.add_argument(
        "--oracle-binder-ch",
        "--oracle_binder_ch",
        dest="oracle_binder_ch",
        default="B",
        help="Binder chain in oracle prediction.pdb files (default: B).",
    )
    parser.add_argument(
        "--oracle-target-ch",
        "--oracle_target_ch",
        dest="oracle_target_ch",
        default="A",
        help="Target chain in oracle prediction.pdb files (default: A).",
    )
    parser.add_argument(
        "--relaxed-binder-ch",
        "--relaxed_binder_ch",
        dest="relaxed_binder_ch",
        default="A",
        help="Binder chain in relaxed.pdb files (default: A).",
    )
    parser.add_argument(
        "--relaxed-target-ch",
        "--relaxed_target_ch",
        dest="relaxed_target_ch",
        default="B",
        help="Target chain in relaxed.pdb files (default: B).",
    )
    parser.add_argument(
        "--round",
        type=int,
        default=4,
        help="Round to analyze (default: 4).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure resolution (default: 300).",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=None,
        help=(
            "Number of bins to use for both histograms. "
            "Defaults to automatic bin selection."
        ),
    )
    parser.add_argument(
        "--hotspot-residues",
        type=parse_hotspot_residues,
        default=None,
        help=(
            "Oracle target-chain residue or inclusive residue range used for "
            "hotspot-distance classification, for example 23 or 20-25."
        ),
    )
    args = parser.parse_args()

    for argument_name in (
        "oracle_binder_ch",
        "oracle_target_ch",
        "relaxed_binder_ch",
        "relaxed_target_ch",
    ):
        value = getattr(args, argument_name)
        if len(value) != 1:
            parser.error(f"--{argument_name.replace('_', '-')} must be one character")
    if args.round < 1:
        parser.error("--round must be at least 1")
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")
    if args.bins is not None and args.bins < 1:
        parser.error("--bins must be at least 1")
    return args


def natural_key(value: str) -> list[object]:
    return [int(piece) if piece.isdigit() else piece.lower() for piece in re.split(r"(\d+)", value)]


def find_oracle_designs(oracle_dir: Path, round_number: int) -> dict[str, Path]:
    if not oracle_dir.is_dir():
        raise ValueError(f"Oracle directory does not exist: {oracle_dir}")

    round_name = f"round_{round_number}"
    designs: dict[str, Path] = {}
    for prediction_path in oracle_dir.glob(f"*/{round_name}/prediction.pdb"):
        design_name = prediction_path.parent.parent.name
        designs[design_name] = prediction_path

    if not designs:
        raise ValueError(
            f"No */{round_name}/prediction.pdb files found below {oracle_dir}"
        )
    return designs


def infer_design_name(pdb_path: Path, design_names: Iterable[str]) -> str:
    design_names = set(design_names)
    for parent in pdb_path.parents:
        if parent.name in design_names:
            return parent.name
    if len(design_names) == 1:
        return next(iter(design_names))
    raise ValueError(
        "Could not infer the design name from the single --relaxed PDB path; "
        "its parent directories must include one of the oracle design names."
    )


def find_relaxed_structures(
    relaxed_input: Path,
    design_names: Iterable[str],
    round_number: int,
) -> dict[str, Path]:
    design_names = set(design_names)
    if relaxed_input.is_file():
        if relaxed_input.suffix.lower() != ".pdb":
            raise ValueError(f"--relaxed is not a PDB file: {relaxed_input}")
        design_name = infer_design_name(relaxed_input, design_names)
        return {design_name: relaxed_input}

    if not relaxed_input.is_dir():
        raise ValueError(f"Relaxed input does not exist: {relaxed_input}")

    round_name = f"round_{round_number}"
    relaxed_paths: dict[str, Path] = {}
    for pdb_path in relaxed_input.glob(f"*/{round_name}/relaxed.pdb"):
        design_name = pdb_path.parent.parent.name
        if design_name in design_names:
            relaxed_paths[design_name] = pdb_path

    # Also support --relaxed pointing at one design or directly at its round.
    direct_candidates = (
        relaxed_input / round_name / "relaxed.pdb",
        relaxed_input / "relaxed.pdb",
    )
    for pdb_path in direct_candidates:
        if pdb_path.is_file():
            design_name = infer_design_name(pdb_path, design_names)
            relaxed_paths[design_name] = pdb_path

    if not relaxed_paths:
        raise ValueError(
            f"No */{round_name}/relaxed.pdb files matching oracle designs found "
            f"below {relaxed_input}"
        )
    return relaxed_paths


def read_ca_atoms(pdb_path: Path, chain_id: str) -> struc.AtomArray:
    pdb_file = pdb.PDBFile.read(str(pdb_path))
    structure = pdb_file.get_structure(model=1)
    ca_atoms = structure[
        (structure.chain_id == chain_id)
        & (structure.atom_name == "CA")
    ]

    if len(ca_atoms) == 0:
        raise ValueError(f"Chain {chain_id!r} has no C-alpha atoms in {pdb_path}")
    return ca_atoms


def calculate_hotspot_distance(
    oracle_pdb: Path,
    target_chain: str,
    binder_chain: str,
    hotspot_residues: tuple[int, ...],
) -> tuple[float, str]:
    pdb_file = pdb.PDBFile.read(str(oracle_pdb))
    structure = pdb_file.get_structure(model=1)
    heavy_atom_mask = structure.element != "H"
    hotspot = structure[
        (structure.chain_id == target_chain)
        & np.isin(structure.res_id, hotspot_residues)
        & heavy_atom_mask
    ]
    binder = structure[
        (structure.chain_id == binder_chain)
        & heavy_atom_mask
    ]

    if len(hotspot) == 0:
        residue_description = ", ".join(str(residue) for residue in hotspot_residues)
        raise ValueError(
            f"No heavy atoms found for hotspot residue(s) {residue_description} "
            f"on oracle target chain {target_chain!r} in {oracle_pdb}"
        )
    if len(binder) == 0:
        raise ValueError(
            f"Oracle binder chain {binder_chain!r} has no heavy atoms in {oracle_pdb}"
        )

    distances = struc.distance(
        binder.coord[:, None, :],
        hotspot.coord[None, :, :],
    )
    minimum_distance = float(np.min(distances))
    classification = HOTSPOT_CATEGORIES[0 if minimum_distance <= 4.0 else 1]
    return minimum_distance, classification


def load_ipae(metrics_path: Path) -> float:
    if not metrics_path.is_file():
        raise ValueError(f"Missing metrics file: {metrics_path}")
    with metrics_path.open(encoding="utf-8") as metrics_file:
        metrics = json.load(metrics_file)
    if "i_pae_normalized" not in metrics:
        raise ValueError(f"Missing i_pae_normalized in {metrics_path}")
    value = float(metrics["i_pae_normalized"])
    if not math.isfinite(value):
        raise ValueError(f"Non-finite i_pae_normalized in {metrics_path}")
    return value


def analyze_design(
    design_name: str,
    oracle_pdb: Path,
    relaxed_pdb: Path,
    ipae: float | str,
    args: argparse.Namespace,
) -> dict[str, object]:
    oracle_target = read_ca_atoms(oracle_pdb, args.oracle_target_ch)
    oracle_binder = read_ca_atoms(oracle_pdb, args.oracle_binder_ch)
    relaxed_target = read_ca_atoms(relaxed_pdb, args.relaxed_target_ch)
    relaxed_binder = read_ca_atoms(relaxed_pdb, args.relaxed_binder_ch)

    if len(relaxed_target) != len(oracle_target):
        raise ValueError(
            "Target C-alpha count differs: "
            f"relaxed={len(relaxed_target)}, oracle={len(oracle_target)}"
        )
    if len(relaxed_binder) != len(oracle_binder):
        raise ValueError(
            "Binder C-alpha count differs: "
            f"relaxed={len(relaxed_binder)}, oracle={len(oracle_binder)}"
        )

    aligned_target, transformation = struc.superimpose(
        oracle_target,
        relaxed_target,
    )
    aligned_binder = transformation.apply(relaxed_binder)
    metrics_path = oracle_pdb.with_name("metrics.json")
    if args.hotspot_residues is None:
        minimum_hotspot_distance: float | str = ""
        hotspot_classification = ""
    else:
        minimum_hotspot_distance, hotspot_classification = (
            calculate_hotspot_distance(
                oracle_pdb,
                args.oracle_target_ch,
                args.oracle_binder_ch,
                args.hotspot_residues,
            )
        )

    return {
        "design": design_name,
        "status": "ok",
        "binder_ca_rmsd_angstrom": float(struc.rmsd(oracle_binder, aligned_binder)),
        "target_alignment_ca_rmsd_angstrom": float(
            struc.rmsd(oracle_target, aligned_target)
        ),
        "i_pae_normalized": ipae,
        "binder_ca_count": len(oracle_binder),
        "target_ca_count": len(oracle_target),
        "binder_residue_name_mismatches": sum(
            first != second
            for first, second in zip(
                relaxed_binder.res_name, oracle_binder.res_name
            )
        ),
        "target_residue_name_mismatches": sum(
            first != second
            for first, second in zip(
                relaxed_target.res_name, oracle_target.res_name
            )
        ),
        "minimum_hotspot_distance_angstrom": minimum_hotspot_distance,
        "hotspot_classification": hotspot_classification,
        "relaxed_pdb": str(relaxed_pdb.resolve()),
        "oracle_pdb": str(oracle_pdb.resolve()),
        "metrics_json": str(metrics_path.resolve()),
        "error": "",
    }


def summary_statistics(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {
            "count": 0,
            "mean": None,
            "standard_deviation": None,
            "median": None,
            "minimum": None,
            "q1": None,
            "q3": None,
            "maximum": None,
        }
    array = np.asarray(values, dtype=float)
    return {
        "count": int(array.size),
        "mean": float(np.mean(array)),
        "standard_deviation": float(np.std(array, ddof=1)) if array.size > 1 else 0.0,
        "median": float(np.median(array)),
        "minimum": float(np.min(array)),
        "q1": float(np.percentile(array, 25)),
        "q3": float(np.percentile(array, 75)),
        "maximum": float(np.max(array)),
    }


def histogram_bins(sample_count: int, requested_bins: int | None) -> str | int:
    if requested_bins is not None:
        return requested_bins
    return "auto" if sample_count > 1 else 1


def plot_histogram(
    values: list[float],
    xlabel: str,
    title: str,
    output_path: Path,
    color: str,
    dpi: int,
    bins: int | None,
) -> None:
    if not values:
        return
    figure, axis = plt.subplots(figsize=(7, 5))
    axis.hist(
        values,
        bins=histogram_bins(len(values), bins),
        color=color,
        edgecolor="white",
        linewidth=0.8,
    )
    mean = float(np.mean(values))
    median = float(np.median(values))
    axis.axvline(mean, color="#b2182b", linestyle="--", linewidth=1.5, label=f"Mean: {mean:.3f}")
    axis.axvline(median, color="#2166ac", linestyle=":", linewidth=1.8, label=f"Median: {median:.3f}")
    axis.set(xlabel=xlabel, ylabel="Design count", title=title)
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def plot_scatter(
    rows: list[dict[str, object]],
    output_path: Path,
    interactive_output_path: Path,
    dpi: int,
) -> None:
    paired_rows = [
        row
        for row in rows
        if row["status"] == "ok"
        and row["binder_ca_rmsd_angstrom"] != ""
        and row["i_pae_normalized"] != ""
    ]
    if not paired_rows:
        return
    x_values = [float(row["binder_ca_rmsd_angstrom"]) for row in paired_rows]
    y_values = [float(row["i_pae_normalized"]) for row in paired_rows]
    hotspot_enabled = paired_rows[0]["hotspot_classification"] != ""

    figure, axis = plt.subplots(figsize=(6, 5))
    if hotspot_enabled:
        for category in HOTSPOT_CATEGORIES:
            category_rows = [
                row
                for row in paired_rows
                if row["hotspot_classification"] == category
            ]
            axis.scatter(
                [float(row["binder_ca_rmsd_angstrom"]) for row in category_rows],
                [float(row["i_pae_normalized"]) for row in category_rows],
                s=34,
                alpha=0.75,
                color=HOTSPOT_COLORS[category],
                edgecolors="white",
                linewidths=0.4,
                label=category,
            )
        axis.legend(title="Hotspot proximity", frameon=False)
    else:
        axis.scatter(x_values, y_values, s=34, alpha=0.75, color="#7b3294", edgecolors="white", linewidths=0.4)
    axis.set(
        xlabel="Target-aligned binder Cα RMSD (Å)",
        ylabel="Normalized iPAE",
        title="Oracle confidence versus structural agreement",
    )
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)

    designs = [str(row["design"]) for row in paired_rows]
    if hotspot_enabled:
        interactive_figure = go.Figure()
        for category in HOTSPOT_CATEGORIES:
            category_rows = [
                row
                for row in paired_rows
                if row["hotspot_classification"] == category
            ]
            interactive_figure.add_trace(
                go.Scatter(
                    x=[
                        float(row["binder_ca_rmsd_angstrom"])
                        for row in category_rows
                    ],
                    y=[float(row["i_pae_normalized"]) for row in category_rows],
                    mode="markers",
                    name=category,
                    text=[str(row["design"]) for row in category_rows],
                    customdata=[
                        [float(row["minimum_hotspot_distance_angstrom"])]
                        for row in category_rows
                    ],
                    marker={"color": HOTSPOT_COLORS[category]},
                    hovertemplate=(
                        "<b>%{text}</b><br>"
                        "Target-aligned binder Cα RMSD (Å): %{x}<br>"
                        "Normalized iPAE: %{y}<br>"
                        "Minimum hotspot distance (Å): %{customdata[0]:.3f}"
                        "<extra>%{fullData.name}</extra>"
                    ),
                )
            )
    else:
        interactive_figure = go.Figure(
            data=go.Scatter(
                x=x_values,
                y=y_values,
                mode="markers",
                text=designs,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "Target-aligned binder Cα RMSD (Å): %{x}<br>"
                    "Normalized iPAE: %{y}<extra></extra>"
                ),
            )
        )
    interactive_figure.update_layout(
        title="Oracle confidence versus structural agreement",
        xaxis_title="Target-aligned binder Cα RMSD (Å)",
        yaxis_title="Normalized iPAE",
    )
    if hotspot_enabled:
        interactive_figure.update_layout(legend_title_text="Hotspot proximity")
    interactive_figure.write_html(
        interactive_output_path,
        include_plotlyjs=True,
    )


def write_results(rows: list[dict[str, object]], output_path: Path) -> None:
    fieldnames = [
        "design",
        "status",
        "binder_ca_rmsd_angstrom",
        "target_alignment_ca_rmsd_angstrom",
        "i_pae_normalized",
        "binder_ca_count",
        "target_ca_count",
        "binder_residue_name_mismatches",
        "target_residue_name_mismatches",
        "minimum_hotspot_distance_angstrom",
        "hotspot_classification",
        "relaxed_pdb",
        "oracle_pdb",
        "metrics_json",
        "error",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_arguments()
    try:
        oracle_designs = find_oracle_designs(args.oracle, args.round)
        relaxed_structures = find_relaxed_structures(
            args.relaxed, oracle_designs, args.round
        )
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    design_names = (
        relaxed_structures.keys() if args.relaxed.is_file() else oracle_designs.keys()
    )
    for design_name in sorted(design_names, key=natural_key):
        oracle_pdb = oracle_designs[design_name]
        metrics_path = oracle_pdb.with_name("metrics.json")
        try:
            ipae: float | str = load_ipae(metrics_path)
            metrics_error = ""
        except (OSError, ValueError, json.JSONDecodeError) as error:
            ipae = ""
            metrics_error = str(error)
        relaxed_pdb = relaxed_structures.get(design_name)
        if relaxed_pdb is None:
            rows.append(
                {
                    "design": design_name,
                    "status": "error",
                    "binder_ca_rmsd_angstrom": "",
                    "target_alignment_ca_rmsd_angstrom": "",
                    "i_pae_normalized": ipae,
                    "binder_ca_count": "",
                    "target_ca_count": "",
                    "binder_residue_name_mismatches": "",
                    "target_residue_name_mismatches": "",
                    "minimum_hotspot_distance_angstrom": "",
                    "hotspot_classification": "",
                    "relaxed_pdb": "",
                    "oracle_pdb": str(oracle_pdb.resolve()),
                    "metrics_json": str(metrics_path.resolve()),
                    "error": "; ".join(
                        value
                        for value in (
                            f"No matching round_{args.round}/relaxed.pdb",
                            metrics_error,
                        )
                        if value
                    ),
                }
            )
            continue
        try:
            row = analyze_design(design_name, oracle_pdb, relaxed_pdb, ipae, args)
            if metrics_error:
                row["error"] = metrics_error
            rows.append(row)
        except (OSError, ValueError, json.JSONDecodeError, np.linalg.LinAlgError) as error:
            rows.append(
                {
                    "design": design_name,
                    "status": "error",
                    "binder_ca_rmsd_angstrom": "",
                    "target_alignment_ca_rmsd_angstrom": "",
                    "i_pae_normalized": ipae,
                    "binder_ca_count": "",
                    "target_ca_count": "",
                    "binder_residue_name_mismatches": "",
                    "target_residue_name_mismatches": "",
                    "minimum_hotspot_distance_angstrom": "",
                    "hotspot_classification": "",
                    "relaxed_pdb": str(relaxed_pdb.resolve()),
                    "oracle_pdb": str(oracle_pdb.resolve()),
                    "metrics_json": str(metrics_path.resolve()),
                    "error": "; ".join(
                        value for value in (str(error), metrics_error) if value
                    ),
                }
            )
            print(f"Warning: {design_name}: {error}", file=sys.stderr)

    successful_rows = [row for row in rows if row["status"] == "ok"]
    rmsd_values = [float(row["binder_ca_rmsd_angstrom"]) for row in successful_rows]
    ipae_values = [
        float(row["i_pae_normalized"])
        for row in rows
        if row["i_pae_normalized"] != ""
    ]

    write_results(rows, args.output_dir / "results.csv")
    summary = {
        "analysis": "binder C-alpha RMSD after target C-alpha superposition",
        "round": args.round,
        "chains": {
            "oracle_binder": args.oracle_binder_ch,
            "oracle_target": args.oracle_target_ch,
            "relaxed_binder": args.relaxed_binder_ch,
            "relaxed_target": args.relaxed_target_ch,
        },
        "designs_discovered": len(rows),
        "designs_analyzed": len(successful_rows),
        "designs_failed": len(rows) - len(successful_rows),
        "binder_ca_rmsd_angstrom": summary_statistics(rmsd_values),
        "i_pae_normalized": summary_statistics(ipae_values),
    }
    with (args.output_dir / "summary.json").open("w", encoding="utf-8") as output_file:
        json.dump(summary, output_file, indent=2)
        output_file.write("\n")

    plot_histogram(
        rmsd_values,
        "Target-aligned binder Cα RMSD (Å)",
        "Binder structural agreement",
        args.output_dir / "rmsd_histogram.png",
        "#4c78a8",
        args.dpi,
        args.bins,
    )
    plot_histogram(
        ipae_values,
        "Normalized iPAE",
        "Oracle interface confidence",
        args.output_dir / "ipae_histogram.png",
        "#f58518",
        args.dpi,
        args.bins,
    )
    plot_scatter(
        rows,
        args.output_dir / "rmsd_vs_ipae.png",
        args.output_dir / "rmsd_vs_ipae.html",
        args.dpi,
    )

    print(f"Analyzed {len(successful_rows)} of {len(rows)} designs.")
    if rmsd_values:
        print(f"Mean binder C-alpha RMSD: {np.mean(rmsd_values):.3f} A")
    if ipae_values:
        print(f"Mean normalized iPAE: {np.mean(ipae_values):.3f}")
    print(f"Wrote results to {args.output_dir.resolve()}")
    return 0 if successful_rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
