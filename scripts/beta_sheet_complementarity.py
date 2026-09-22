#!/usr/bin/env python3

"""
Analyze inter- and intramolecular beta-sheet complementarity in RFdiffusion designs.

The script uses Biotite to read and validate PDB structures and mkdssp to assign
beta-sheet topology.

A design is counted as having beta-sheet complementarity when the configured
minimum number of DSSP E-state residue pairs belong to the same beta ladder.

Outputs:
    beta_sheet_summary.csv
    inter_beta_pairs.csv
    intra_beta_pairs.csv
    dssp/*.dssp
    beta_sheet_complementarity.png
"""

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import biotite.structure.io.pdb as pdb


@dataclass
class DsspResidue:
    dssp_index: int
    pdb_residue_number: int
    insertion_code: str
    chain_id: str
    amino_acid: str
    secondary_structure: str
    bridge_label_1: str
    bridge_label_2: str
    bridge_partner_1: int
    bridge_partner_2: int


@dataclass
class BetaPair:
    residue_1: DsspResidue
    residue_2: DsspResidue
    ladder_label: str
    orientation: str


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Analyze inter- and intramolecular beta-sheet complementarity "
            "for a directory of PDB structures."
        )
    )

    parser.add_argument(
        "--input-pdbs",
        required=True,
        type=Path,
        help="Directory containing PDB files.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for analysis outputs.",
    )

    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=None,
        help="Directory for figures. Defaults to --output-dir.",
    )

    parser.add_argument(
        "--target-chain",
        required=True,
        help="PDB chain ID for the target.",
    )

    parser.add_argument(
        "--target-interface",
        nargs="+",
        default=None,
        help=(
            "Target residue numbers. Accepts individual numbers, comma-separated "
            "numbers, and ranges such as 45-52. If omitted, the full target chain "
            "is used."
        ),
    )

    parser.add_argument(
        "--binder-chain",
        required=True,
        help="PDB chain ID for the binder.",
    )

    parser.add_argument(
        "--min-inter-beta-pairs",
        type=int,
        default=3,
        help=(
            "Minimum number of DSSP beta-paired residue pairs in the same "
            "ladder for inter beta-sheet complementarity (default: 3)."
        ),
    )

    parser.add_argument(
        "--min-intra-beta-pairs",
        type=int,
        default=2,
        help=(
            "Minimum number of DSSP beta-paired residue pairs in the same "
            "ladder for intra beta-sheet complementarity (default: 2)."
        ),
    )

    return parser.parse_args()


def parse_target_interface(values):
    if values is None:
        return None

    residue_numbers = set()

    for value in values:
        pieces = value.split(",")

        for piece in pieces:
            piece = piece.strip()

            if piece == "":
                continue

            if "-" in piece:
                range_parts = piece.split("-")

                if len(range_parts) != 2:
                    raise ValueError(
                        f"Could not parse target-interface value: {piece}"
                    )

                start = int(range_parts[0])
                end = int(range_parts[1])

                if end < start:
                    raise ValueError(
                        f"Target-interface range ends before it starts: {piece}"
                    )

                for residue_number in range(start, end + 1):
                    residue_numbers.add(residue_number)

            else:
                residue_numbers.add(int(piece))

    if len(residue_numbers) == 0:
        raise ValueError("No target-interface residues were parsed.")

    return residue_numbers


def residue_label(residue_number, insertion_code):
    if insertion_code:
        return f"{residue_number}{insertion_code}"

    return str(residue_number)


def get_chain_residues(pdb_path, chain_id):
    pdb_file = pdb.PDBFile.read(str(pdb_path))
    structure = pdb_file.get_structure(model=1)

    chain_mask = structure.chain_id == chain_id
    protein_mask = ~structure.hetero
    ca_mask = structure.atom_name == "CA"

    chain_ca_atoms = structure[chain_mask & protein_mask & ca_mask]

    if len(chain_ca_atoms) == 0:
        raise ValueError(
            f"Chain {chain_id} was not found or contains no protein CA atoms."
        )

    residues = []

    for res_id, ins_code in zip(
        chain_ca_atoms.res_id,
        chain_ca_atoms.ins_code,
    ):
        residue = (int(res_id), str(ins_code).strip())

        if residue not in residues:
            residues.append(residue)

    return residues


def write_dssp_input_pdb(input_pdb_path, output_pdb_path):
    # Biotite fills missing element fields from atom names when the PDB is read.
    input_pdb = pdb.PDBFile.read(str(input_pdb_path))
    structure = input_pdb.get_structure(model=1)

    output_pdb = pdb.PDBFile()
    output_pdb.set_structure(structure)
    output_pdb.write(str(output_pdb_path))


def run_dssp(pdb_path, dssp_path):
    # RFdiffusion PDB files may have blank element columns. mkdssp rejects
    # these files, so a temporary PDB with populated element fields is used.
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_pdb_path = Path(temporary_directory) / pdb_path.name

        write_dssp_input_pdb(
            pdb_path,
            temporary_pdb_path,
        )

        command = [
            "mkdssp",
            "--output-format=dssp",
            str(temporary_pdb_path),
            str(dssp_path),
        ]

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

    if result.returncode != 0:
        error_message = result.stderr.strip()

        if error_message == "":
            error_message = "mkdssp returned a non-zero exit code."

        raise RuntimeError(error_message)

    if not dssp_path.exists():
        raise RuntimeError("mkdssp did not create the expected DSSP output file.")


def parse_integer_field(text):
    text = text.strip()

    if text == "":
        return 0

    return int(text)


def parse_dssp_file(dssp_path):
    records = []
    residue_section_started = False

    with dssp_path.open() as handle:
        for line in handle:
            if "RESIDUE AA STRUCTURE" in line:
                residue_section_started = True
                continue

            if not residue_section_started:
                continue

            if len(line) < 34:
                continue

            # DSSP uses an exclamation mark for chain-break records.
            if len(line) > 13 and line[13] == "!":
                continue

            try:
                dssp_index = int(line[0:5])
                pdb_residue_number = int(line[5:10])
            except ValueError:
                continue

            insertion_code = line[10].strip()
            chain_id = line[11].strip()
            amino_acid = line[13].strip()

            secondary_structure = line[16].strip()

            if secondary_structure == "":
                secondary_structure = "-"

            bridge_label_1 = line[23].strip()
            bridge_label_2 = line[24].strip()

            bridge_partner_1 = parse_integer_field(line[25:29])
            bridge_partner_2 = parse_integer_field(line[29:33])

            record = DsspResidue(
                dssp_index=dssp_index,
                pdb_residue_number=pdb_residue_number,
                insertion_code=insertion_code,
                chain_id=chain_id,
                amino_acid=amino_acid,
                secondary_structure=secondary_structure,
                bridge_label_1=bridge_label_1,
                bridge_label_2=bridge_label_2,
                bridge_partner_1=bridge_partner_1,
                bridge_partner_2=bridge_partner_2,
            )

            records.append(record)

    if len(records) == 0:
        raise RuntimeError(
            f"No residue records were parsed from DSSP output: {dssp_path}"
        )

    return records


def bridge_orientation(bridge_label):
    if bridge_label == "":
        return "unknown"

    if bridge_label.islower():
        return "parallel"

    if bridge_label.isupper():
        return "antiparallel"

    return "unknown"


def extract_beta_pairs(records):
    records_by_index = {}

    for record in records:
        records_by_index[record.dssp_index] = record

    beta_pairs = []
    seen_pairs = set()

    for record in records:
        # E means that the residue participates in a beta ladder.
        if record.secondary_structure != "E":
            continue

        partner_data = [
            (record.bridge_partner_1, record.bridge_label_1),
            (record.bridge_partner_2, record.bridge_label_2),
        ]

        for partner_index, bridge_label in partner_data:
            if partner_index == 0:
                continue

            if bridge_label == "":
                continue

            partner = records_by_index.get(partner_index)

            if partner is None:
                continue

            if partner.secondary_structure != "E":
                continue

            pair_key = tuple(
                sorted([record.dssp_index, partner.dssp_index])
            )

            if pair_key in seen_pairs:
                continue

            seen_pairs.add(pair_key)

            beta_pair = BetaPair(
                residue_1=record,
                residue_2=partner,
                ladder_label=bridge_label,
                orientation=bridge_orientation(bridge_label),
            )

            beta_pairs.append(beta_pair)

    return beta_pairs


def split_inter_and_intra_pairs(
    beta_pairs,
    target_chain,
    binder_chain,
    target_interface_numbers,
):
    inter_pairs = []
    intra_pairs = []

    for pair in beta_pairs:
        chain_1 = pair.residue_1.chain_id
        chain_2 = pair.residue_2.chain_id

        is_target_binder_pair = (
            chain_1 == target_chain and chain_2 == binder_chain
        ) or (
            chain_1 == binder_chain and chain_2 == target_chain
        )

        if is_target_binder_pair:
            if chain_1 == target_chain:
                target_residue = pair.residue_1
                binder_residue = pair.residue_2
            else:
                target_residue = pair.residue_2
                binder_residue = pair.residue_1

            if (
                target_residue.pdb_residue_number
                in target_interface_numbers
            ):
                inter_pairs.append(
                    (
                        target_residue,
                        binder_residue,
                        pair.ladder_label,
                        pair.orientation,
                    )
                )

        if chain_1 == binder_chain and chain_2 == binder_chain:
            intra_pairs.append(pair)

    return inter_pairs, intra_pairs


def has_beta_ladder(pair_records, minimum_pair_count):
    pairs_by_ladder = {}

    for pair_record in pair_records:
        if isinstance(pair_record, BetaPair):
            ladder_label = pair_record.ladder_label
        else:
            ladder_label = pair_record[2]

        if ladder_label not in pairs_by_ladder:
            pairs_by_ladder[ladder_label] = 0

        pairs_by_ladder[ladder_label] += 1

    for pair_count in pairs_by_ladder.values():
        if pair_count >= minimum_pair_count:
            return True

    return False


def write_summary_csv(summary_rows, output_path):
    fieldnames = [
        "design",
        "status",
        "inter_beta_complementarity",
        "inter_target_beta_pair_count",
        "intra_beta_complementarity",
        "intra_binder_beta_pair_count",
        "error",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in summary_rows:
            writer.writerow(row)


def write_inter_pair_csv(rows, output_path):
    fieldnames = [
        "design",
        "target_chain",
        "target_residue",
        "binder_chain",
        "binder_residue",
        "ladder",
        "orientation",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def write_intra_pair_csv(rows, output_path):
    fieldnames = [
        "design",
        "binder_chain",
        "binder_residue_1",
        "binder_residue_2",
        "ladder",
        "orientation",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def choose_tick_positions(number_of_items):
    if number_of_items <= 30:
        return np.arange(number_of_items)

    step = int(np.ceil(number_of_items / 30))
    return np.arange(0, number_of_items, step)


def make_figure(
    heatmap_percent,
    target_labels,
    binder_labels,
    inter_positive_count,
    intra_positive_count,
    analyzed_count,
    figure_path,
):
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(14, 7),
        constrained_layout=True,
    )

    heatmap_axis = axes[0]

    heatmap_for_plot = np.ma.masked_equal(
        heatmap_percent,
        0,
    )

    cmap = plt.colormaps["Blues"].copy()
    cmap.set_bad("white")

    image = heatmap_axis.imshow(
        heatmap_for_plot,
        cmap=cmap,
        aspect="auto",
        origin="lower",
        vmin=0,
        vmax=100,
    )

    heatmap_axis.set_title("Target-binder beta pairing")
    heatmap_axis.set_xlabel("Binder residue")
    heatmap_axis.set_ylabel("Target residue")

    binder_tick_positions = choose_tick_positions(len(binder_labels))
    target_tick_positions = choose_tick_positions(len(target_labels))

    heatmap_axis.set_xticks(binder_tick_positions)
    heatmap_axis.set_xticklabels(
        [binder_labels[index] for index in binder_tick_positions],
        rotation=90,
    )

    heatmap_axis.set_yticks(target_tick_positions)
    heatmap_axis.set_yticklabels(
        [target_labels[index] for index in target_tick_positions]
    )

    binder_cell_boundaries = np.arange(
        -0.5,
        len(binder_labels),
        1,
    )
    target_cell_boundaries = np.arange(
        -0.5,
        len(target_labels),
        1,
    )

    heatmap_axis.set_xticks(
        binder_cell_boundaries,
        minor=True,
    )
    heatmap_axis.set_yticks(
        target_cell_boundaries,
        minor=True,
    )
    heatmap_axis.grid(
        which="minor",
        color="lightgray",
        linewidth=0.5,
    )
    heatmap_axis.tick_params(
        which="minor",
        bottom=False,
        left=False,
    )

    colorbar = figure.colorbar(
        image,
        ax=heatmap_axis,
        fraction=0.046,
        pad=0.04,
    )
    colorbar.set_label("Designs with beta pair (%)")

    summary_axis = axes[1]

    if analyzed_count == 0:
        inter_percent = 0.0
        intra_percent = 0.0
    else:
        inter_percent = 100.0 * inter_positive_count / analyzed_count
        intra_percent = 100.0 * intra_positive_count / analyzed_count

    category_names = [
        "Inter\n(target interface)",
        "Intra\n(binder)",
    ]
    percentages = [
        inter_percent,
        intra_percent,
    ]
    counts = [
        inter_positive_count,
        intra_positive_count,
    ]

    bars = summary_axis.bar(
        category_names,
        percentages,
    )

    summary_axis.set_ylim(0, 100)
    summary_axis.set_ylabel("Designs with beta-sheet complementarity (%)")
    summary_axis.set_title(
        f"Binder beta-sheet complementarity\nN = {analyzed_count}"
    )

    for bar, percentage, count in zip(
        bars,
        percentages,
        counts,
    ):
        summary_axis.text(
            bar.get_x() + bar.get_width() / 2,
            percentage + 2,
            f"{percentage:.1f}%\n({count}/{analyzed_count})",
            ha="center",
            va="bottom",
        )

    figure.savefig(
        figure_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(figure)


def main():
    args = parse_arguments()

    print(
        "Minimum inter beta-pair threshold: "
        f"{args.min_inter_beta_pairs}"
    )
    print(
        "Minimum intra beta-pair threshold: "
        f"{args.min_intra_beta_pairs}"
    )

    if args.target_chain == args.binder_chain:
        print(
            "ERROR: --target-chain and --binder-chain must be different.",
            file=sys.stderr,
        )
        return 1

    if shutil.which("mkdssp") is None:
        print(
            "ERROR: mkdssp was not found in the active environment.",
            file=sys.stderr,
        )
        return 1

    if not args.input_pdbs.is_dir():
        print(
            f"ERROR: Input directory does not exist: {args.input_pdbs}",
            file=sys.stderr,
        )
        return 1

    pdb_files = sorted(args.input_pdbs.glob("*.pdb"))

    if len(pdb_files) == 0:
        print(
            f"ERROR: No PDB files found in: {args.input_pdbs}",
            file=sys.stderr,
        )
        return 1

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if args.figure_dir is None:
        figure_dir = args.output_dir
    else:
        figure_dir = args.figure_dir
        figure_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

    dssp_dir = args.output_dir / "dssp"
    dssp_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:
        requested_target_interface = parse_target_interface(
            args.target_interface
        )
    except ValueError as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

    first_pdb = pdb_files[0]

    try:
        target_residues = get_chain_residues(
            first_pdb,
            args.target_chain,
        )
        binder_residues = get_chain_residues(
            first_pdb,
            args.binder_chain,
        )
    except Exception as error:
        print(
            f"ERROR: Could not read reference chains from {first_pdb.name}: {error}",
            file=sys.stderr,
        )
        return 1

    target_residue_numbers = set()

    for residue_number, insertion_code in target_residues:
        target_residue_numbers.add(residue_number)

    if requested_target_interface is None:
        target_interface_numbers = target_residue_numbers
        selected_target_residues = target_residues
    else:
        missing_residues = (
            requested_target_interface - target_residue_numbers
        )

        if len(missing_residues) > 0:
            missing_text = ", ".join(
                str(number)
                for number in sorted(missing_residues)
            )

            print(
                (
                    "ERROR: The following target-interface residues were not "
                    f"found in chain {args.target_chain}: {missing_text}"
                ),
                file=sys.stderr,
            )
            return 1

        target_interface_numbers = requested_target_interface
        selected_target_residues = []

        for residue in target_residues:
            if residue[0] in target_interface_numbers:
                selected_target_residues.append(residue)

    target_labels = []

    for residue_number, insertion_code in selected_target_residues:
        target_labels.append(
            residue_label(
                residue_number,
                insertion_code,
            )
        )

    binder_labels = []

    for residue_number, insertion_code in binder_residues:
        binder_labels.append(
            residue_label(
                residue_number,
                insertion_code,
            )
        )

    target_index = {}
    binder_index = {}

    for index, residue in enumerate(selected_target_residues):
        target_index[residue] = index

    for index, residue in enumerate(binder_residues):
        binder_index[residue] = index

    heatmap_counts = np.zeros(
        (
            len(selected_target_residues),
            len(binder_residues),
        ),
        dtype=int,
    )

    summary_rows = []
    inter_pair_rows = []
    intra_pair_rows = []

    analyzed_count = 0
    inter_positive_count = 0
    intra_positive_count = 0

    for pdb_path in pdb_files:
        design_name = pdb_path.stem

        print(f"Analyzing {pdb_path.name}")

        summary_row = {
            "design": design_name,
            "status": "ok",
            "inter_beta_complementarity": "",
            "inter_target_beta_pair_count": "",
            "intra_beta_complementarity": "",
            "intra_binder_beta_pair_count": "",
            "error": "",
        }

        try:
            current_target_residues = get_chain_residues(
                pdb_path,
                args.target_chain,
            )
            current_binder_residues = get_chain_residues(
                pdb_path,
                args.binder_chain,
            )

            if current_binder_residues != binder_residues:
                raise ValueError(
                    "Binder residue numbering differs from the first PDB file."
                )

            current_target_numbers = set()

            for residue_number, insertion_code in current_target_residues:
                current_target_numbers.add(residue_number)

            missing_target_residues = (
                target_interface_numbers - current_target_numbers
            )

            if len(missing_target_residues) > 0:
                raise ValueError(
                    "One or more target-interface residues are missing."
                )

            dssp_path = dssp_dir / f"{design_name}.dssp"

            run_dssp(
                pdb_path,
                dssp_path,
            )

            dssp_records = parse_dssp_file(dssp_path)

            beta_pairs = extract_beta_pairs(dssp_records)

            inter_pairs, intra_pairs = split_inter_and_intra_pairs(
                beta_pairs,
                args.target_chain,
                args.binder_chain,
                target_interface_numbers,
            )

            inter_positive = has_beta_ladder(
                inter_pairs,
                args.min_inter_beta_pairs,
            )
            intra_positive = has_beta_ladder(
                intra_pairs,
                args.min_intra_beta_pairs,
            )

            analyzed_count += 1

            if inter_positive:
                inter_positive_count += 1

            if intra_positive:
                intra_positive_count += 1

            summary_row["inter_beta_complementarity"] = int(
                inter_positive
            )
            summary_row["inter_target_beta_pair_count"] = len(
                inter_pairs
            )
            summary_row["intra_beta_complementarity"] = int(
                intra_positive
            )
            summary_row["intra_binder_beta_pair_count"] = len(
                intra_pairs
            )

            pairs_seen_in_this_design = set()

            for (
                target_residue,
                binder_residue,
                ladder_label,
                orientation,
            ) in inter_pairs:
                target_key = (
                    target_residue.pdb_residue_number,
                    target_residue.insertion_code,
                )
                binder_key = (
                    binder_residue.pdb_residue_number,
                    binder_residue.insertion_code,
                )

                matrix_pair_key = (
                    target_key,
                    binder_key,
                )

                if matrix_pair_key not in pairs_seen_in_this_design:
                    target_position = target_index.get(target_key)
                    binder_position = binder_index.get(binder_key)

                    if (
                        target_position is not None
                        and binder_position is not None
                    ):
                        heatmap_counts[
                            target_position,
                            binder_position,
                        ] += 1

                    pairs_seen_in_this_design.add(matrix_pair_key)

                inter_pair_rows.append(
                    {
                        "design": design_name,
                        "target_chain": args.target_chain,
                        "target_residue": residue_label(
                            target_residue.pdb_residue_number,
                            target_residue.insertion_code,
                        ),
                        "binder_chain": args.binder_chain,
                        "binder_residue": residue_label(
                            binder_residue.pdb_residue_number,
                            binder_residue.insertion_code,
                        ),
                        "ladder": ladder_label,
                        "orientation": orientation,
                    }
                )

            for pair in intra_pairs:
                intra_pair_rows.append(
                    {
                        "design": design_name,
                        "binder_chain": args.binder_chain,
                        "binder_residue_1": residue_label(
                            pair.residue_1.pdb_residue_number,
                            pair.residue_1.insertion_code,
                        ),
                        "binder_residue_2": residue_label(
                            pair.residue_2.pdb_residue_number,
                            pair.residue_2.insertion_code,
                        ),
                        "ladder": pair.ladder_label,
                        "orientation": pair.orientation,
                    }
                )

        except Exception as error:
            summary_row["status"] = "failed"
            summary_row["error"] = str(error)

            print(
                f"WARNING: {pdb_path.name} failed: {error}",
                file=sys.stderr,
            )

        summary_rows.append(summary_row)

    if analyzed_count == 0:
        print(
            "ERROR: No PDB files were successfully analyzed.",
            file=sys.stderr,
        )

        write_summary_csv(
            summary_rows,
            args.output_dir / "beta_sheet_summary.csv",
        )

        return 1

    heatmap_percent = (
        heatmap_counts.astype(float)
        / analyzed_count
        * 100.0
    )

    write_summary_csv(
        summary_rows,
        args.output_dir / "beta_sheet_summary.csv",
    )

    write_inter_pair_csv(
        inter_pair_rows,
        args.output_dir / "inter_beta_pairs.csv",
    )

    write_intra_pair_csv(
        intra_pair_rows,
        args.output_dir / "intra_beta_pairs.csv",
    )

    heatmap_csv_path = (
        args.output_dir
        / "inter_beta_pair_frequency_percent.csv"
    )

    with heatmap_csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow(
            ["target_residue"] + binder_labels
        )

        for row_index, target_label in enumerate(target_labels):
            row = [target_label]

            for value in heatmap_percent[row_index]:
                row.append(f"{value:.3f}")

            writer.writerow(row)

    figure_path = (
        figure_dir
        / "beta_sheet_complementarity.png"
    )

    make_figure(
        heatmap_percent,
        target_labels,
        binder_labels,
        inter_positive_count,
        intra_positive_count,
        analyzed_count,
        figure_path,
    )

    inter_percent = (
        100.0
        * inter_positive_count
        / analyzed_count
    )

    intra_percent = (
        100.0
        * intra_positive_count
        / analyzed_count
    )

    print()
    print(f"Successfully analyzed: {analyzed_count}")
    print(
        "Inter beta-sheet complementarity: "
        f"{inter_percent:.1f}% "
        f"({inter_positive_count}/{analyzed_count})"
    )
    print(
        "Intra beta-sheet complementarity: "
        f"{intra_percent:.1f}% "
        f"({intra_positive_count}/{analyzed_count})"
    )
    print(f"Figure: {figure_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
