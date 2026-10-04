#!/usr/bin/env python3

"""
Calculate cyclic backbone torsion angles for RFD3 macrocycles.

For every residue, this script records:

    phi   = C(previous) - N - CA - C
    psi   = N - CA - C - N(next)
    omega = CA - C - N(next) - CA(next)

The first and last residues are connected explicitly. This is important for
cyclic peptides: a linear-chain helper would leave the first phi and final
psi/omega undefined.

Outputs:

    ramachandran/all_structures_torsions.csv
    ramachandran/validation_summary.json

Run this script with the existing Biotite conda environment.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import statistics
from dataclasses import dataclass
from pathlib import Path

import biotite
import biotite.structure as struc
import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure


# =============================================================================
# Paths and filename convention
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
    / "ramachandran"
)

INPUT_FILENAME_PATTERN = re.compile(
    r"macrocycle_(?P<length>\d+)"
    r"_(?P<batch>\d+)"
    r"_model_(?P<model>\d+)"
    r"\.cif\.gz$"
)

BACKBONE_ATOM_NAMES = ("N", "CA", "C")


# =============================================================================
# Data structures
# =============================================================================


@dataclass(frozen=True)
class DesignFile:
    """One generated macrocycle structure."""

    path: Path
    peptide_length: int
    batch: int
    model: int

    @property
    def generation_order(self) -> tuple[int, int]:
        return self.batch, self.model


@dataclass(frozen=True)
class TorsionRecord:
    """Backbone measurements for one residue."""

    peptide_length: int
    batch: int
    model: int
    filename: str
    chain_id: str
    residue_index: int
    residue_id: int
    residue_name: str
    next_residue_name: str
    ramachandran_class: str
    peptide_bond_class: str
    phi_deg: float
    psi_deg: float
    omega_deg: float
    c_to_next_n_angstrom: float
    is_closure_bond: bool


@dataclass
class LengthValidation:
    """Measurements collected while processing one peptide length."""

    peptide_length: int
    structure_count: int
    residue_count: int
    closure_distances: list[float]
    peptide_bond_distances: list[float]
    phi_angles: list[float]
    psi_angles: list[float]
    omega_angles: list[float]
    cis_proline_count: int
    cis_general_count: int
    twisted_proline_count: int
    twisted_general_count: int


# =============================================================================
# Arguments
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate cyclic phi, psi, and omega angles with Biotite."
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
        help=f"Directory for torsion results. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=[10, 12],
        help="Peptide lengths to analyze. Default: 10 12",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing torsion and validation files.",
    )

    return parser.parse_args()


# =============================================================================
# Input discovery
# =============================================================================


def discover_designs(
    input_dir: Path,
    requested_lengths: set[int],
) -> list[DesignFile]:
    """Find and identify generated macrocycle structures."""

    designs: list[DesignFile] = []

    for path in input_dir.glob("*.cif.gz"):
        match = INPUT_FILENAME_PATTERN.search(path.name)

        if match is None:
            continue

        peptide_length = int(match.group("length"))

        if peptide_length not in requested_lengths:
            continue

        designs.append(
            DesignFile(
                path=path.resolve(),
                peptide_length=peptide_length,
                batch=int(match.group("batch")),
                model=int(match.group("model")),
            )
        )

    designs.sort(
        key=lambda design: (
            design.peptide_length,
            design.generation_order,
        )
    )

    seen: set[tuple[int, int, int]] = set()

    for design in designs:
        identifier = (
            design.peptide_length,
            design.batch,
            design.model,
        )

        if identifier in seen:
            raise ValueError(
                "Duplicate structure identifier: "
                f"length={design.peptide_length}, "
                f"batch={design.batch}, model={design.model}"
            )

        seen.add(identifier)

    return designs


# =============================================================================
# Structure parsing and validation
# =============================================================================


def read_structure(path: Path) -> struc.AtomArray:
    """Read one gzip-compressed mmCIF model with Biotite."""

    with gzip.open(path, "rt") as handle:
        cif_file = CIFFile.read(handle)

    return get_structure(
        cif_file,
        model=1,
        altloc="first",
    )


def get_backbone_coordinates(
    residue: struc.AtomArray,
    *,
    path: Path,
    residue_index: int,
) -> dict[str, np.ndarray]:
    """Return the unique N, CA, and C coordinates for one residue."""

    coordinates: dict[str, np.ndarray] = {}

    for atom_name in BACKBONE_ATOM_NAMES:
        matching_coordinates = residue.coord[
            residue.atom_name == atom_name
        ]

        if len(matching_coordinates) != 1:
            raise ValueError(
                f"Expected exactly one {atom_name} atom in residue "
                f"{residue_index} of {path.name}; found "
                f"{len(matching_coordinates)}."
            )

        coordinates[atom_name] = matching_coordinates[0]

    return coordinates


def classify_ramachandran_residue(
    residue_name: str,
    next_residue_name: str,
    incoming_omega_degrees: float,
) -> str:
    """Assign one of the six residue classes used by MolProbity."""

    if residue_name == "GLY":
        return "glycine"

    if residue_name == "PRO":
        if classify_omega(incoming_omega_degrees) == "cis":
            return "cis-proline"

        return "trans-proline"

    if next_residue_name == "PRO":
        return "pre-proline"

    if residue_name in {"ILE", "VAL"}:
        return "isoleucine or valine"

    return "general"


def classify_omega(omega_degrees: float) -> str:
    """Classify an omega angle using the MolProbity 30-degree convention."""

    if abs(omega_degrees) <= 30.0:
        return "cis"

    if abs(abs(omega_degrees) - 180.0) <= 30.0:
        return "trans"

    return "twisted"


def calculate_torsions(
    design: DesignFile,
) -> list[TorsionRecord]:
    """Calculate all cyclic backbone torsions for one structure."""

    atoms = read_structure(design.path)

    chain_ids = sorted({str(chain_id) for chain_id in atoms.chain_id})

    if len(chain_ids) != 1:
        raise ValueError(
            f"Expected one chain in {design.path.name}; found {chain_ids}."
        )

    chain_id = chain_ids[0]
    amino_acid_atoms = atoms[struc.filter_amino_acids(atoms)]
    residues = list(struc.residue_iter(amino_acid_atoms))

    if len(residues) != design.peptide_length:
        raise ValueError(
            f"Filename identifies {design.peptide_length} residues but "
            f"Biotite found {len(residues)} in {design.path.name}."
        )

    backbone_coordinates = [
        get_backbone_coordinates(
            residue,
            path=design.path,
            residue_index=index,
        )
        for index, residue in enumerate(residues, start=1)
    ]

    records: list[TorsionRecord] = []

    for zero_based_index, residue in enumerate(residues):
        # Allows cyclic wrapping for index 0
        # e.g., (0 - 1) % 10 = 9
        previous_index = (zero_based_index - 1) % len(residues)
        # Allow syclic wrapping for index 10
        # e.g., (10 + 1) % 10 = 1
        next_index = (zero_based_index + 1) % len(residues)

        previous_atoms = backbone_coordinates[previous_index]
        current_atoms = backbone_coordinates[zero_based_index]
        next_atoms = backbone_coordinates[next_index]

        phi_radians = struc.dihedral(
            previous_atoms["C"],
            current_atoms["N"],
            current_atoms["CA"],
            current_atoms["C"],
        )

        psi_radians = struc.dihedral(
            current_atoms["N"],
            current_atoms["CA"],
            current_atoms["C"],
            next_atoms["N"],
        )

        omega_radians = struc.dihedral(
            current_atoms["CA"],
            current_atoms["C"],
            next_atoms["N"],
            next_atoms["CA"],
        )

        phi_degrees = float(np.degrees(phi_radians))
        psi_degrees = float(np.degrees(psi_radians))
        omega_degrees = float(np.degrees(omega_radians))

        incoming_omega_radians = struc.dihedral(
            previous_atoms["CA"],
            previous_atoms["C"],
            current_atoms["N"],
            current_atoms["CA"],
        )

        incoming_omega_degrees = float(
            np.degrees(incoming_omega_radians)
        )

        measured_angles = (
            phi_degrees,
            psi_degrees,
            omega_degrees,
        )

        if not all(math.isfinite(angle) for angle in measured_angles):
            raise ValueError(
                f"Non-finite backbone angle in residue "
                f"{zero_based_index + 1} of {design.path.name}."
            )

        residue_name = str(residue.res_name[0])
        next_residue_name = str(residues[next_index].res_name[0])
        peptide_bond_class = (
            "X-Pro"
            if next_residue_name == "PRO"
            else "non-Pro"
        )

        c_to_next_n = float(
            struc.distance(
                current_atoms["C"],
                next_atoms["N"],
            )
        )

        records.append(
            TorsionRecord(
                peptide_length=design.peptide_length,
                batch=design.batch,
                model=design.model,
                filename=design.path.name,
                chain_id=chain_id,
                residue_index=zero_based_index + 1,
                residue_id=int(residue.res_id[0]),
                residue_name=residue_name,
                next_residue_name=next_residue_name,
                ramachandran_class=classify_ramachandran_residue(
                    residue_name,
                    next_residue_name,
                    incoming_omega_degrees,
                ),
                peptide_bond_class=peptide_bond_class,
                phi_deg=phi_degrees,
                psi_deg=psi_degrees,
                omega_deg=omega_degrees,
                c_to_next_n_angstrom=c_to_next_n,
                is_closure_bond=(next_index == 0),
            )
        )

    return records


# =============================================================================
# Output and summary helpers
# =============================================================================


TORSION_FIELDNAMES = [
    "peptide_length",
    "batch",
    "model",
    "filename",
    "chain_id",
    "residue_index",
    "residue_id",
    "residue_name",
    "next_residue_name",
    "ramachandran_class",
    "peptide_bond_class",
    "phi_deg",
    "psi_deg",
    "omega_deg",
    "c_to_next_n_angstrom",
    "is_closure_bond",
]


def make_length_validation(peptide_length: int) -> LengthValidation:
    """Create an empty validation accumulator."""

    return LengthValidation(
        peptide_length=peptide_length,
        structure_count=0,
        residue_count=0,
        closure_distances=[],
        peptide_bond_distances=[],
        phi_angles=[],
        psi_angles=[],
        omega_angles=[],
        cis_proline_count=0,
        cis_general_count=0,
        twisted_proline_count=0,
        twisted_general_count=0,
    )


def update_validation(
    validation: LengthValidation,
    records: list[TorsionRecord],
) -> None:
    """Add one structure's measurements to the validation accumulator."""

    validation.structure_count += 1
    validation.residue_count += len(records)

    for record in records:
        validation.phi_angles.append(record.phi_deg)
        validation.psi_angles.append(record.psi_deg)
        validation.omega_angles.append(record.omega_deg)
        validation.peptide_bond_distances.append(
            record.c_to_next_n_angstrom
        )

        if record.is_closure_bond:
            validation.closure_distances.append(
                record.c_to_next_n_angstrom
            )

        omega_class = classify_omega(record.omega_deg)
        is_x_pro = record.peptide_bond_class == "X-Pro"

        if omega_class == "cis":
            if is_x_pro:
                validation.cis_proline_count += 1
            else:
                validation.cis_general_count += 1

        if omega_class == "twisted":
            if is_x_pro:
                validation.twisted_proline_count += 1
            else:
                validation.twisted_general_count += 1


def describe_values(values: list[float]) -> dict[str, float]:
    """Return simple descriptive statistics for a numeric list."""

    return {
        "minimum": min(values),
        "percentile_0.1": float(np.percentile(values, 0.1)),
        "percentile_1": float(np.percentile(values, 1.0)),
        "maximum": max(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "percentile_99": float(np.percentile(values, 99.0)),
        "percentile_99.9": float(np.percentile(values, 99.9)),
    }


def broad_peptide_bond_distance_screen(
    values: list[float],
) -> dict[str, float | int | str]:
    """Count conspicuous C-N distances using a deliberately broad interval."""

    lower_bound = 1.2
    upper_bound = 1.5
    outside_count = sum(
        value < lower_bound or value > upper_bound
        for value in values
    )

    return {
        "purpose": (
            "simple diagnostic screen, not a restraint-dictionary z-score"
        ),
        "lower_bound_angstrom": lower_bound,
        "upper_bound_angstrom": upper_bound,
        "outside_count": outside_count,
        "outside_percent": 100.0 * outside_count / len(values),
    }


def validation_to_dict(
    validation: LengthValidation,
) -> dict[str, object]:
    """Convert one length's validation accumulator to JSON-ready data."""

    expected_residue_count = (
        validation.structure_count
        * validation.peptide_length
    )
    trans_deviations = [
        abs(abs(angle) - 180.0)
        for angle in validation.omega_angles
        if classify_omega(angle) == "trans"
    ]

    return {
        "peptide_length": validation.peptide_length,
        "structure_count": validation.structure_count,
        "residue_count": validation.residue_count,
        "expected_residue_count": expected_residue_count,
        "all_expected_residues_present": (
            validation.residue_count == expected_residue_count
        ),
        "all_torsions_finite": True,
        "closure_c_to_n_angstrom": describe_values(
            validation.closure_distances
        ),
        "closure_c_to_n_broad_screen": broad_peptide_bond_distance_screen(
            validation.closure_distances
        ),
        "all_peptide_bond_c_to_n_angstrom": describe_values(
            validation.peptide_bond_distances
        ),
        "all_peptide_bond_c_to_n_broad_screen": (
            broad_peptide_bond_distance_screen(
                validation.peptide_bond_distances
            )
        ),
        "phi_degrees": describe_values(validation.phi_angles),
        "psi_degrees": describe_values(validation.psi_angles),
        "omega_degrees": describe_values(validation.omega_angles),
        "trans_omega_deviation_from_180_degrees": describe_values(
            trans_deviations
        ),
        "omega_geometry": {
            "definition": {
                "cis": "abs(omega) <= 30 degrees",
                "trans": "abs(abs(omega) - 180) <= 30 degrees",
                "twisted": "all remaining omega angles",
            },
            "cis_x_pro": validation.cis_proline_count,
            "cis_non_pro": validation.cis_general_count,
            "twisted_x_pro": validation.twisted_proline_count,
            "twisted_non_pro": validation.twisted_general_count,
        },
    }


# =============================================================================
# Main analysis
# =============================================================================


def main() -> None:
    args = parse_arguments()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist:\n{input_dir}"
        )

    requested_lengths = set(args.lengths)

    if not requested_lengths:
        raise ValueError("At least one peptide length is required.")

    designs = discover_designs(
        input_dir,
        requested_lengths,
    )

    if not designs:
        raise FileNotFoundError(
            f"No matching .cif.gz structures found in:\n{input_dir}"
        )

    found_lengths = {design.peptide_length for design in designs}
    missing_lengths = requested_lengths - found_lengths

    if missing_lengths:
        raise FileNotFoundError(
            "No structures were found for peptide length(s): "
            + ", ".join(str(length) for length in sorted(missing_lengths))
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    torsion_path = output_dir / "all_structures_torsions.csv"
    summary_path = output_dir / "validation_summary.json"

    existing_outputs = [
        path
        for path in (torsion_path, summary_path)
        if path.exists()
    ]

    if existing_outputs and not args.overwrite:
        raise FileExistsError(
            "Output already exists. Use --overwrite to replace it:\n"
            + "\n".join(str(path) for path in existing_outputs)
        )

    validations = {
        peptide_length: make_length_validation(peptide_length)
        for peptide_length in sorted(requested_lengths)
    }

    temporary_torsion_path = output_dir / ".all_structures_torsions.csv.tmp"

    print(f"Input directory : {input_dir}")
    print(f"Structures      : {len(designs):,}")
    print(
        "Lengths         : "
        + ", ".join(str(length) for length in sorted(requested_lengths))
    )
    print(f"Biotite version : {biotite.__version__}")
    print(f"Output directory: {output_dir}")

    try:
        with temporary_torsion_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=TORSION_FIELDNAMES,
                lineterminator="\n",
            )
            writer.writeheader()

            for index, design in enumerate(designs, start=1):
                records = calculate_torsions(design)

                for record in records:
                    writer.writerow(
                        {
                            field_name: getattr(record, field_name)
                            for field_name in TORSION_FIELDNAMES
                        }
                    )

                update_validation(
                    validations[design.peptide_length],
                    records,
                )

                if index % 1000 == 0 or index == len(designs):
                    print(
                        f"Processed {index:,}/{len(designs):,} structures",
                        flush=True,
                    )

        temporary_torsion_path.replace(torsion_path)

    except Exception:
        if temporary_torsion_path.exists():
            temporary_torsion_path.unlink()
        raise

    summary = {
        "analysis": "cyclic macrocycle backbone torsions",
        "input_directory": str(input_dir),
        "torsion_table": str(torsion_path),
        "software": {
            "biotite": biotite.__version__,
            "numpy": np.__version__,
        },
        "method": {
            "coordinate_parser": "biotite.structure.io.pdbx.CIFFile",
            "dihedral_function": "biotite.structure.dihedral",
            "angle_unit": "degrees",
            "cyclic_wraparound": True,
            "phi_atoms": "C(previous)-N-CA-C",
            "psi_atoms": "N-CA-C-N(next)",
            "omega_atoms": "CA-C-N(next)-CA(next)",
            "ramachandran_residue_classes": [
                "general",
                "glycine",
                "cis-proline",
                "trans-proline",
                "pre-proline",
                "isoleucine or valine",
            ],
        },
        "results": [
            validation_to_dict(validations[length])
            for length in sorted(validations)
        ],
    }

    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(f"Saved: {torsion_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()
