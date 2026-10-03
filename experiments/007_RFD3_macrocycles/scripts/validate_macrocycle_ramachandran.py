#!/usr/bin/env python3

"""
Compare macrocycle phi/psi angles with MolProbity Top8000 distributions.

This script reads the torsion table made by
``calculate_macrocycle_torsions.py`` and evaluates every residue against the
same six reference distributions and cutoffs used by ``phenix.ramalyze``.

The bundled reference files are unmodified copies of the official Richardson
Lab ``rotarama_data`` repository. The interpolation follows cctbx's
``NDimTable.valueAt()`` behavior: linear interpolation between 2-degree bin
centers with periodic wrapping at -180/180 degrees.

Outputs:

    ramachandran/ramachandran_validation.csv
    ramachandran/ramachandran_quality.csv
    ramachandran/ramachandran_validation_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy
from scipy.ndimage import map_coordinates


# =============================================================================
# Paths and MolProbity definitions
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
DEFAULT_OUTPUT_DIR = DEFAULT_RAMACHANDRAN_DIR
DEFAULT_REFERENCE_DIR = (
    EXPERIMENT_DIR
    / "reference_data"
    / "molprobity_top8000"
)

REFERENCE_FILENAMES = {
    "general": "rama8000-general-noGPIVpreP.data",
    "glycine": "rama8000-gly-sym.data",
    "cis-proline": "rama8000-cispro.data",
    "trans-proline": "rama8000-transpro.data",
    "pre-proline": "rama8000-prepro-noGP.data",
    "isoleucine or valine": "rama8000-ileval-nopreP.data",
}

RESIDUE_CLASS_ORDER = list(REFERENCE_FILENAMES)

FAVORED_CUTOFF = 0.02

OUTLIER_CUTOFFS = {
    "general": 0.0005,
    "cis-proline": 0.0020,
    "glycine": 0.0010,
    "trans-proline": 0.0010,
    "pre-proline": 0.0010,
    "isoleucine or valine": 0.0010,
}

QUALITY_ORDER = ["favored", "allowed", "outlier"]


# =============================================================================
# Data structures
# =============================================================================


@dataclass
class ValidationRecord:
    """Identifiers, torsions, and MolProbity result for one residue."""

    peptide_length: int
    batch: int
    model: int
    filename: str
    residue_index: int
    residue_id: int
    residue_name: str
    ramachandran_class: str
    phi_deg: float
    psi_deg: float
    top8000_score: float = 0.0
    evaluation: str = ""

    @property
    def structure_key(self) -> tuple[int, int, int]:
        return self.peptide_length, self.batch, self.model


# =============================================================================
# Arguments
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate macrocycle phi/psi angles against Top8000."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help=f"Torsion CSV. Default: {DEFAULT_INPUT_PATH}",
    )

    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help=f"Directory containing rama8000 data. Default: {DEFAULT_REFERENCE_DIR}",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for validation results. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing validation outputs.",
    )

    return parser.parse_args()


# =============================================================================
# Input loading
# =============================================================================


def read_torsions(input_path: Path) -> list[ValidationRecord]:
    """Read the columns needed for Top8000 evaluation."""

    records: list[ValidationRecord] = []

    with input_path.open() as handle:
        reader = csv.DictReader(handle)

        required_columns = {
            "peptide_length",
            "batch",
            "model",
            "filename",
            "residue_index",
            "residue_id",
            "residue_name",
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
            residue_class = row["ramachandran_class"]

            if residue_class not in REFERENCE_FILENAMES:
                raise ValueError(
                    f"Unknown Ramachandran class {residue_class!r} in "
                    f"{input_path}."
                )

            records.append(
                ValidationRecord(
                    peptide_length=int(row["peptide_length"]),
                    batch=int(row["batch"]),
                    model=int(row["model"]),
                    filename=row["filename"],
                    residue_index=int(row["residue_index"]),
                    residue_id=int(row["residue_id"]),
                    residue_name=row["residue_name"],
                    ramachandran_class=residue_class,
                    phi_deg=float(row["phi_deg"]),
                    psi_deg=float(row["psi_deg"]),
                )
            )

    return records


def load_top8000_grid(data_path: Path) -> np.ndarray:
    """Load one sparse cctbx NDimTable text file into a 180 x 180 grid."""

    grid = np.zeros((180, 180), dtype=float)

    with data_path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.split()

            if len(fields) != 3:
                raise ValueError(
                    f"Unexpected row in {data_path.name}: {line.rstrip()}"
                )

            phi, psi, value = (float(field) for field in fields)

            phi_index = int((phi + 179.0) / 2.0)
            psi_index = int((psi + 179.0) / 2.0)

            if not (0 <= phi_index < 180 and 0 <= psi_index < 180):
                raise ValueError(
                    f"Grid coordinate outside -179..179 in {data_path.name}: "
                    f"phi={phi}, psi={psi}"
                )

            grid[phi_index, psi_index] = value

    return grid


# =============================================================================
# MolProbity evaluation
# =============================================================================


def interpolate_scores(
    grid: np.ndarray,
    phi: np.ndarray,
    psi: np.ndarray,
) -> np.ndarray:
    """Interpolate a periodic Top8000 grid at arbitrary phi/psi angles."""

    # Grid centers are -179, -177, ..., 179 degrees. Convert each angle to
    # a fractional array index, then let SciPy perform bilinear interpolation
    # with periodic wrapping across both dihedral boundaries.
    phi_coordinates = (phi + 179.0) / 2.0
    psi_coordinates = (psi + 179.0) / 2.0

    return map_coordinates(
        grid,
        [phi_coordinates, psi_coordinates],
        order=1,
        mode="grid-wrap",
        prefilter=False,
    )


def classify_score(residue_class: str, score: float) -> str:
    """Apply the cutoffs from cctbx ``ramalyze.evalScore()``."""

    if score >= FAVORED_CUTOFF:
        return "favored"

    if score >= OUTLIER_CUTOFFS[residue_class]:
        return "allowed"

    return "outlier"


def evaluate_records(
    records: list[ValidationRecord],
    reference_dir: Path,
) -> None:
    """Evaluate all records in place, one reference class at a time."""

    for residue_class, filename in REFERENCE_FILENAMES.items():
        data_path = reference_dir / filename

        if not data_path.is_file():
            raise FileNotFoundError(
                f"Missing Top8000 reference file:\n{data_path}"
            )

        selected_indices = [
            index
            for index, record in enumerate(records)
            if record.ramachandran_class == residue_class
        ]

        if not selected_indices:
            continue

        phi = np.asarray(
            [records[index].phi_deg for index in selected_indices],
            dtype=float,
        )
        psi = np.asarray(
            [records[index].psi_deg for index in selected_indices],
            dtype=float,
        )

        scores = interpolate_scores(
            load_top8000_grid(data_path),
            phi,
            psi,
        )

        for index, score in zip(selected_indices, scores):
            record = records[index]
            record.top8000_score = float(score)
            record.evaluation = classify_score(
                residue_class,
                record.top8000_score,
            )


# =============================================================================
# Outputs
# =============================================================================


def write_detailed_results(
    records: list[ValidationRecord],
    output_path: Path,
) -> None:
    """Write one Top8000 score and classification per residue."""

    fieldnames = [
        "peptide_length",
        "batch",
        "model",
        "filename",
        "residue_index",
        "residue_id",
        "residue_name",
        "ramachandran_class",
        "phi_deg",
        "psi_deg",
        "top8000_score",
        "evaluation",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()

        for record in records:
            writer.writerow(
                {
                    field_name: getattr(record, field_name)
                    for field_name in fieldnames
                }
            )


def count_quality(
    records: list[ValidationRecord],
) -> dict[tuple[int, str], Counter[str]]:
    """Count favored, allowed, and outlier residues by length and class."""

    counts: dict[tuple[int, str], Counter[str]] = defaultdict(Counter)

    for record in records:
        counts[(record.peptide_length, "all")][record.evaluation] += 1
        counts[
            (record.peptide_length, record.ramachandran_class)
        ][record.evaluation] += 1

    return counts


def write_quality_table(
    records: list[ValidationRecord],
    output_path: Path,
) -> None:
    """Write aggregate residue-level validation counts and percentages."""

    counts = count_quality(records)
    peptide_lengths = sorted({record.peptide_length for record in records})

    with output_path.open("w", newline="") as handle:
        fieldnames = [
            "peptide_length",
            "ramachandran_class",
            "n_residues",
            "n_favored",
            "percent_favored",
            "n_allowed",
            "percent_allowed",
            "n_outlier",
            "percent_outlier",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()

        for peptide_length in peptide_lengths:
            for residue_class in ["all", *RESIDUE_CLASS_ORDER]:
                quality_counts = counts[(peptide_length, residue_class)]
                total = sum(quality_counts.values())

                if total == 0:
                    continue

                row: dict[str, int | float | str] = {
                    "peptide_length": peptide_length,
                    "ramachandran_class": residue_class,
                    "n_residues": total,
                }

                for quality in QUALITY_ORDER:
                    count = quality_counts[quality]
                    row[f"n_{quality}"] = count
                    row[f"percent_{quality}"] = 100.0 * count / total

                writer.writerow(row)


def make_summary(
    records: list[ValidationRecord],
    *,
    input_path: Path,
    reference_dir: Path,
    detailed_path: Path,
    quality_path: Path,
) -> dict[str, object]:
    """Create a machine-readable validation and provenance summary."""

    counts = count_quality(records)
    peptide_lengths = sorted({record.peptide_length for record in records})
    structures_by_length: dict[int, set[tuple[int, int, int]]] = defaultdict(set)
    outliers_by_structure: Counter[tuple[int, int, int]] = Counter()

    for record in records:
        structures_by_length[record.peptide_length].add(record.structure_key)

        if record.evaluation == "outlier":
            outliers_by_structure[record.structure_key] += 1

    results = []

    for peptide_length in peptide_lengths:
        quality_counts = counts[(peptide_length, "all")]
        residue_count = sum(quality_counts.values())
        structures = structures_by_length[peptide_length]
        structures_with_outliers = sum(
            outliers_by_structure[key] > 0
            for key in structures
        )

        results.append(
            {
                "peptide_length": peptide_length,
                "structure_count": len(structures),
                "residue_count": residue_count,
                "favored": {
                    "count": quality_counts["favored"],
                    "percent": 100.0
                    * quality_counts["favored"]
                    / residue_count,
                },
                "allowed": {
                    "count": quality_counts["allowed"],
                    "percent": 100.0
                    * quality_counts["allowed"]
                    / residue_count,
                },
                "outlier": {
                    "count": quality_counts["outlier"],
                    "percent": 100.0
                    * quality_counts["outlier"]
                    / residue_count,
                },
                "structures_with_at_least_one_outlier": {
                    "count": structures_with_outliers,
                    "percent": 100.0
                    * structures_with_outliers
                    / len(structures),
                },
            }
        )

    return {
        "analysis": "MolProbity Top8000 Ramachandran validation",
        "input_torsion_table": str(input_path),
        "detailed_output": str(detailed_path),
        "quality_table": str(quality_path),
        "software": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "reference": {
            "name": "Richardson Lab Top8000 Ramachandran distributions",
            "directory": str(reference_dir),
            "source": "https://github.com/rlabduke/rotarama_data",
            "commit": "76aae74c6e1f834775f8df4700c79602c6a8b9ec",
            "license": "CC BY 4.0",
            "residue_classes": RESIDUE_CLASS_ORDER,
        },
        "method": {
            "interpolation": (
                "bilinear interpolation between 2-degree bin centers with "
                "periodic wrapping, matching cctbx NDimTable.valueAt"
            ),
            "favored_cutoff": FAVORED_CUTOFF,
            "outlier_cutoffs": OUTLIER_CUTOFFS,
            "cutoff_source": "cctbx mmtbx.validation.ramalyze.evalScore",
            "cyclic_wraparound_torsions_included": True,
        },
        "results": results,
    }


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    args = parse_arguments()

    input_path = args.input.resolve()
    reference_dir = args.reference_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"Torsion CSV does not exist:\n{input_path}")

    if not reference_dir.is_dir():
        raise NotADirectoryError(
            f"Top8000 reference directory does not exist:\n{reference_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    detailed_path = output_dir / "ramachandran_validation.csv"
    quality_path = output_dir / "ramachandran_quality.csv"
    summary_path = output_dir / "ramachandran_validation_summary.json"
    output_paths = [detailed_path, quality_path, summary_path]

    existing_paths = [path for path in output_paths if path.exists()]

    if existing_paths and not args.overwrite:
        raise FileExistsError(
            "Output already exists. Use --overwrite to replace it:\n"
            + "\n".join(str(path) for path in existing_paths)
        )

    print(f"Input          : {input_path}")
    print(f"Top8000 data   : {reference_dir}")

    records = read_torsions(input_path)
    print(f"Residues       : {len(records):,}")

    evaluate_records(records, reference_dir)

    write_detailed_results(records, detailed_path)
    write_quality_table(records, quality_path)

    summary = make_summary(
        records,
        input_path=input_path,
        reference_dir=reference_dir,
        detailed_path=detailed_path,
        quality_path=quality_path,
    )

    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    for path in output_paths:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
