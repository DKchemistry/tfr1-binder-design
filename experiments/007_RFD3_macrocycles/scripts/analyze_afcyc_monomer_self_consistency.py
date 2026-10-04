#!/usr/bin/env python3

"""Analyze the RFD3 -> LigandMPNN -> AfCyc monomer experiment.

The primary endpoint follows RFpeptides Figure 1c:

* evaluate 200 backbones at each peptide length;
* test eight LigandMPNN sequences for each backbone;
* call a backbone successful if at least one sequence has pLDDT > 0.8
  and backbone RMSD < 2.0 Angstrom.

RMSD is calculated over the index-matched backbone-heavy atoms N, CA, C, and O
after least-squares superposition with Biotite.  The lowest RMSD over cyclic
residue shifts is also recorded as a diagnostic, but it does not determine
success: RFpeptides describes cyclic permutations for structural clustering,
not for this self-consistency calculation.

This script reproduces the raw ``all`` statistic from RFpeptides Figure 1c.
The paper's ``unique @ TM 0.5`` statistic additionally requires MaxCluster
clustering across all 1,200 backbones and is intentionally kept separate.

The script only reads an AfCyc prediction after its ``.complete`` marker has
been written.  It is therefore safe to run while the prediction campaign is
still in progress.  Provisional results use only backbones with all eight
predictions.  The final RFpeptides-style success fraction is written only when
every backbone is either complete or has an invalid LigandMPNN FASTA file.
"""

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import biotite.structure as struc
import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure


PEPTIDE_LENGTHS = (8, 10, 12, 14, 16, 18)
BACKBONE_ATOMS = ("N", "CA", "C", "O")
EXPECTED_BACKBONES = 200
EXPECTED_SEQUENCES = 8

SEQUENCE_DIR_PATTERN = re.compile(r"sequence_(\d+)$")

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1] / "monomer_self_consistency"


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--experiment-dir",
        type=Path,
        default=EXPERIMENT_DIR,
        help=f"Monomer experiment directory. Default: {EXPERIMENT_DIR}",
    )
    parser.add_argument(
        "--run-name",
        default="run_01",
        help="LigandMPNN and AfCyc run to analyze. Default: run_01",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: <experiment-dir>/analysis/<run-name>",
    )
    parser.add_argument(
        "--plddt-cutoff",
        type=float,
        default=0.8,
        help="A sequence must be above this pLDDT. Default: 0.8",
    )
    parser.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=2.0,
        help="A sequence must be below this backbone RMSD in Angstrom. Default: 2.0",
    )
    return parser.parse_args()


def model_number(path):
    """Return the integer model number at the end of an RFD3 filename."""

    match = re.search(r"_model_(\d+)\.cif\.gz$", path.name)
    if match is None:
        raise ValueError(f"Could not read model number from {path.name}")
    return int(match.group(1))


def read_cif(path):
    """Read the first model from a gzip-compressed CIF file."""

    with gzip.open(path, mode="rt") as handle:
        return get_structure(CIFFile.read(handle), model=1)


def read_pdb(path):
    """Read the first model from a PDB file."""

    return PDBFile.read(path).get_structure(model=1)


def backbone_coordinates(atom_array, expected_length, source_path):
    """Return one N/CA/C/O coordinate set for each residue in file order."""

    residue_starts = struc.get_residue_starts(
        atom_array,
        add_exclusive_stop=True,
    )
    coordinates = []

    for start, stop in zip(residue_starts[:-1], residue_starts[1:]):
        residue = atom_array[start:stop]
        residue_coordinates = []

        for atom_name in BACKBONE_ATOMS:
            matches = np.where(residue.atom_name == atom_name)[0]

            if len(matches) != 1:
                residue_label = (
                    f"{residue.chain_id[0]}{residue.res_id[0]} "
                    f"{residue.res_name[0]}"
                )
                raise ValueError(
                    f"Expected one {atom_name} atom in {residue_label} of "
                    f"{source_path}; found {len(matches)}"
                )

            residue_coordinates.append(residue.coord[matches[0]])

        coordinates.append(residue_coordinates)

    coordinates = np.asarray(coordinates)

    if len(coordinates) != expected_length:
        raise ValueError(
            f"Expected {expected_length} residues in {source_path}; "
            f"found {len(coordinates)}"
        )

    return coordinates


def superimposed_rmsd(fixed_coordinates, mobile_coordinates):
    """Superimpose two coordinate arrays and return their RMSD."""

    fitted_coordinates, _ = struc.superimpose(
        fixed_coordinates,
        mobile_coordinates,
    )
    return float(struc.rmsd(fixed_coordinates, fitted_coordinates))


def cyclic_backbone_rmsd(source_coordinates, prediction_coordinates):
    """Return the lowest N/CA/C/O RMSD over all cyclic residue shifts."""

    source_flat = source_coordinates.reshape(-1, 3)
    best_rmsd = None
    best_shift = None

    for shift in range(len(prediction_coordinates)):
        shifted_prediction = np.roll(
            prediction_coordinates,
            shift=-shift,
            axis=0,
        )
        prediction_flat = shifted_prediction.reshape(-1, 3)
        rmsd = superimposed_rmsd(source_flat, prediction_flat)

        if best_rmsd is None or rmsd < best_rmsd:
            best_rmsd = rmsd
            best_shift = shift

    return best_rmsd, best_shift


def closure_distance(backbone):
    """Return the distance from the final carbonyl C to the first N."""

    n_index = BACKBONE_ATOMS.index("N")
    c_index = BACKBONE_ATOMS.index("C")
    return float(np.linalg.norm(backbone[-1, c_index] - backbone[0, n_index]))


def read_fasta_records(path):
    """Return FASTA headers and sequences without requiring another package."""

    records = []
    header = None
    sequence_lines = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_lines)))
                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)

    if header is not None:
        records.append((header, "".join(sequence_lines)))

    return records


def ligandmpnn_fasta_is_valid(path, peptide_length):
    """Check that a FASTA contains eight numbered sequences of the right length."""

    if not path.is_file():
        return False, "missing LigandMPNN FASTA"

    sequences = {}

    for header, sequence in read_fasta_records(path):
        match = re.search(r"(?:^|,\s*)id=(\d+)(?:,|$)", header)
        if match is not None:
            sequences[int(match.group(1))] = sequence

    expected_ids = set(range(1, EXPECTED_SEQUENCES + 1))
    if set(sequences) != expected_ids:
        return False, f"sequence IDs are {sorted(sequences)}"

    wrong_lengths = sorted(
        (sequence_id, len(sequence))
        for sequence_id, sequence in sequences.items()
        if len(sequence) != peptide_length
    )
    if wrong_lengths:
        return False, f"wrong sequence lengths: {wrong_lengths}"

    return True, ""


def completed_prediction_dirs(backbone_dir):
    """Return completed sequence directories in sequence-number order."""

    completed = []

    if not backbone_dir.is_dir():
        return completed

    for sequence_dir in backbone_dir.glob("sequence_*"):
        match = SEQUENCE_DIR_PATTERN.fullmatch(sequence_dir.name)
        if match is None or not (sequence_dir / ".complete").is_file():
            continue
        completed.append((int(match.group(1)), sequence_dir))

    return sorted(completed)


def analyze_prediction(
    source_path,
    source_coordinates,
    peptide_length,
    backbone_name,
    sequence_id,
    prediction_dir,
    plddt_cutoff,
    rmsd_cutoff,
):
    """Analyze one completed AfCyc sequence prediction."""

    metrics_path = prediction_dir / "metrics.json"
    prediction_path = prediction_dir / "prediction.pdb"

    if not metrics_path.is_file() or not prediction_path.is_file():
        raise FileNotFoundError(
            f"Completed prediction is missing metrics or structure: {prediction_dir}"
        )

    with metrics_path.open() as handle:
        metrics = json.load(handle)

    if metrics["backbone"] != backbone_name:
        raise ValueError(f"Backbone name mismatch in {metrics_path}")
    if metrics["sequence_id"] != sequence_id:
        raise ValueError(f"Sequence ID mismatch in {metrics_path}")
    if metrics["peptide_length"] != peptide_length:
        raise ValueError(f"Peptide length mismatch in {metrics_path}")

    prediction_coordinates = backbone_coordinates(
        read_pdb(prediction_path),
        expected_length=peptide_length,
        source_path=prediction_path,
    )
    source_flat = source_coordinates.reshape(-1, 3)
    prediction_flat = prediction_coordinates.reshape(-1, 3)
    backbone_rmsd = superimposed_rmsd(source_flat, prediction_flat)

    ca_index = BACKBONE_ATOMS.index("CA")
    ca_rmsd = superimposed_rmsd(
        source_coordinates[:, ca_index],
        prediction_coordinates[:, ca_index],
    )

    cyclic_rmsd, cyclic_shift = cyclic_backbone_rmsd(
        source_coordinates,
        prediction_coordinates,
    )

    plddt = float(metrics["best_plddt"])
    passes_plddt = plddt > plddt_cutoff
    passes_rmsd = backbone_rmsd < rmsd_cutoff

    return {
        "peptide_length": peptide_length,
        "backbone": backbone_name,
        "sequence_id": sequence_id,
        "sequence": metrics["sequence"],
        "afcyc_model": metrics["best_model_by_plddt"],
        "plddt": plddt,
        "backbone_rmsd_angstrom": backbone_rmsd,
        "ca_rmsd_angstrom": ca_rmsd,
        "cyclic_backbone_rmsd_angstrom": cyclic_rmsd,
        "best_cyclic_shift": cyclic_shift,
        "source_closure_distance_angstrom": closure_distance(source_coordinates),
        "prediction_closure_distance_angstrom": closure_distance(
            prediction_coordinates
        ),
        "passes_plddt": passes_plddt,
        "passes_rmsd": passes_rmsd,
        "success": passes_plddt and passes_rmsd,
        "source_path": str(source_path.resolve()),
        "prediction_path": str(prediction_path.resolve()),
    }


def backbone_status(fasta_valid, fasta_reason, prediction_rows):
    """Describe whether a backbone is ready for the final analysis."""

    if not fasta_valid:
        return "invalid_sequence_input", fasta_reason
    if len(prediction_rows) == EXPECTED_SEQUENCES:
        return "complete", ""
    if prediction_rows:
        return "partial", f"{len(prediction_rows)} of 8 predictions complete"
    return "pending", "no predictions complete"


def summarize_backbone(
    peptide_length,
    backbone_name,
    fasta_valid,
    fasta_reason,
    prediction_rows,
):
    """Summarize the eight sequence attempts for one backbone."""

    status, status_detail = backbone_status(
        fasta_valid,
        fasta_reason,
        prediction_rows,
    )
    successful_rows = [row for row in prediction_rows if row["success"]]
    plddt_passes = [row for row in prediction_rows if row["passes_plddt"]]
    rmsd_passes = [row for row in prediction_rows if row["passes_rmsd"]]

    return {
        "peptide_length": peptide_length,
        "backbone": backbone_name,
        "status": status,
        "status_detail": status_detail,
        "completed_sequences": len(prediction_rows),
        "successful_sequences": len(successful_rows),
        "success": bool(successful_rows),
        "first_successful_sequence_id": (
            min(row["sequence_id"] for row in successful_rows)
            if successful_rows
            else ""
        ),
        "max_plddt": (
            max(row["plddt"] for row in prediction_rows)
            if prediction_rows
            else ""
        ),
        "min_backbone_rmsd_angstrom": (
            min(row["backbone_rmsd_angstrom"] for row in prediction_rows)
            if prediction_rows
            else ""
        ),
        "any_sequence_passes_plddt": bool(plddt_passes),
        "any_sequence_passes_rmsd": bool(rmsd_passes),
    }


def summarize_lengths(backbone_rows):
    """Calculate RFpeptides-style and failure-mode summaries by length."""

    rows_by_length = defaultdict(list)
    for row in backbone_rows:
        rows_by_length[row["peptide_length"]].append(row)

    summaries = []

    for peptide_length in PEPTIDE_LENGTHS:
        rows = rows_by_length[peptide_length]
        complete_rows = [row for row in rows if row["status"] == "complete"]
        invalid_rows = [
            row for row in rows if row["status"] == "invalid_sequence_input"
        ]
        successful_rows = [row for row in complete_rows if row["success"]]
        plddt_rows = [
            row for row in complete_rows if row["any_sequence_passes_plddt"]
        ]
        rmsd_rows = [
            row for row in complete_rows if row["any_sequence_passes_rmsd"]
        ]

        campaign_complete = (
            len(complete_rows) + len(invalid_rows) == len(rows)
            and len(rows) == EXPECTED_BACKBONES
        )

        summaries.append(
            {
                "peptide_length": peptide_length,
                "expected_backbones": EXPECTED_BACKBONES,
                "complete_backbones": len(complete_rows),
                "invalid_backbones": len(invalid_rows),
                "pending_or_partial_backbones": (
                    len(rows) - len(complete_rows) - len(invalid_rows)
                ),
                "successful_backbones": len(successful_rows),
                "provisional_success_fraction": (
                    len(successful_rows) / len(complete_rows)
                    if complete_rows
                    else ""
                ),
                "rfpeptides_style_success_fraction": (
                    len(successful_rows) / EXPECTED_BACKBONES
                    if campaign_complete
                    else ""
                ),
                "plddt_only_pass_fraction_complete": (
                    len(plddt_rows) / len(complete_rows)
                    if complete_rows
                    else ""
                ),
                "rmsd_only_pass_fraction_complete": (
                    len(rmsd_rows) / len(complete_rows)
                    if complete_rows
                    else ""
                ),
                "campaign_complete": campaign_complete,
            }
        )

    return summaries


def summarize_attempts(prediction_rows, backbone_rows):
    """Calculate success after the first one through eight sequence attempts."""

    predictions_by_backbone = defaultdict(list)
    for row in prediction_rows:
        key = (row["peptide_length"], row["backbone"])
        predictions_by_backbone[key].append(row)

    summaries = []

    for peptide_length in PEPTIDE_LENGTHS:
        complete_backbones = [
            row
            for row in backbone_rows
            if row["peptide_length"] == peptide_length
            and row["status"] == "complete"
        ]

        for attempt_count in range(1, EXPECTED_SEQUENCES + 1):
            successful_backbones = 0

            for backbone_row in complete_backbones:
                key = (peptide_length, backbone_row["backbone"])
                attempts = predictions_by_backbone[key]
                if any(
                    row["success"] and row["sequence_id"] <= attempt_count
                    for row in attempts
                ):
                    successful_backbones += 1

            summaries.append(
                {
                    "peptide_length": peptide_length,
                    "sequence_attempts": attempt_count,
                    "complete_backbones": len(complete_backbones),
                    "successful_backbones": successful_backbones,
                    "success_fraction": (
                        successful_backbones / len(complete_backbones)
                        if complete_backbones
                        else ""
                    ),
                }
            )

    return summaries


def write_tsv(path, rows):
    """Write dictionaries as a tab-separated table."""

    if not rows:
        raise ValueError(f"No rows available for {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_arguments()
    experiment_dir = args.experiment_dir.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else experiment_dir / "analysis" / args.run_name
    )

    if not 0.0 <= args.plddt_cutoff <= 1.0:
        raise ValueError("--plddt-cutoff must be between zero and one")
    if args.rmsd_cutoff <= 0.0:
        raise ValueError("--rmsd-cutoff must be greater than zero")

    prediction_rows = []
    backbone_rows = []

    for peptide_length in PEPTIDE_LENGTHS:
        length_name = f"length_{peptide_length:02d}"
        source_dir = experiment_dir / "rfd3_backbones" / length_name
        fasta_dir = (
            experiment_dir
            / "ligandmpnn_sequences"
            / args.run_name
            / length_name
            / "seqs"
        )
        prediction_length_dir = (
            experiment_dir
            / "afcyc_predictions"
            / args.run_name
            / length_name
        )

        source_paths = sorted(source_dir.glob("*.cif.gz"), key=model_number)
        if len(source_paths) != EXPECTED_BACKBONES:
            raise ValueError(
                f"Expected {EXPECTED_BACKBONES} backbones in {source_dir}; "
                f"found {len(source_paths)}"
            )

        for source_path in source_paths:
            backbone_name = source_path.name.removesuffix(".cif.gz")
            fasta_path = fasta_dir / f"{backbone_name}.fa"
            fasta_valid, fasta_reason = ligandmpnn_fasta_is_valid(
                fasta_path,
                peptide_length,
            )

            prediction_dirs = completed_prediction_dirs(
                prediction_length_dir / backbone_name
            )
            backbone_predictions = []

            if prediction_dirs:
                source_coordinates = backbone_coordinates(
                    read_cif(source_path),
                    expected_length=peptide_length,
                    source_path=source_path,
                )

                for sequence_id, prediction_dir in prediction_dirs:
                    row = analyze_prediction(
                        source_path=source_path,
                        source_coordinates=source_coordinates,
                        peptide_length=peptide_length,
                        backbone_name=backbone_name,
                        sequence_id=sequence_id,
                        prediction_dir=prediction_dir,
                        plddt_cutoff=args.plddt_cutoff,
                        rmsd_cutoff=args.rmsd_cutoff,
                    )
                    prediction_rows.append(row)
                    backbone_predictions.append(row)

            backbone_rows.append(
                summarize_backbone(
                    peptide_length=peptide_length,
                    backbone_name=backbone_name,
                    fasta_valid=fasta_valid,
                    fasta_reason=fasta_reason,
                    prediction_rows=backbone_predictions,
                )
            )

    prediction_rows.sort(
        key=lambda row: (
            row["peptide_length"],
            model_number(Path(row["source_path"])),
            row["sequence_id"],
        )
    )
    backbone_rows.sort(
        key=lambda row: (
            row["peptide_length"],
            int(row["backbone"].rsplit("_model_", 1)[1]),
        )
    )

    length_rows = summarize_lengths(backbone_rows)
    attempt_rows = summarize_attempts(prediction_rows, backbone_rows)

    write_tsv(output_dir / "prediction_results.tsv", prediction_rows)
    write_tsv(output_dir / "backbone_results.tsv", backbone_rows)
    write_tsv(output_dir / "success_by_length.tsv", length_rows)
    write_tsv(output_dir / "success_by_sequence_attempts.tsv", attempt_rows)

    summary = {
        "method": {
            "reference": "RFpeptides Figure 1c",
            "backbone_atoms": list(BACKBONE_ATOMS),
            "primary_residue_correspondence": "index-matched",
            "cyclic_permutation_rmsd_recorded_as_diagnostic": True,
            "plddt_operator": ">",
            "plddt_cutoff": args.plddt_cutoff,
            "rmsd_operator": "<",
            "rmsd_cutoff_angstrom": args.rmsd_cutoff,
            "expected_backbones_per_length": EXPECTED_BACKBONES,
            "expected_sequences_per_backbone": EXPECTED_SEQUENCES,
        },
        "run_name": args.run_name,
        "completed_predictions_analyzed": len(prediction_rows),
        "length_summaries": length_rows,
    }

    with (output_dir / "analysis_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(f"Analyzed {len(prediction_rows):,} completed predictions")
    print(f"Results: {output_dir}")
    for row in length_rows:
        print(
            f"Length {row['peptide_length']:2d}: "
            f"{row['complete_backbones']:3d} complete backbones, "
            f"{row['successful_backbones']:3d} successful, "
            f"{row['invalid_backbones']} invalid"
        )


if __name__ == "__main__":
    main()
