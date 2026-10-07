#!/usr/bin/env python3
"""Convert RFD3 CIF outputs into PDB inputs for LigandMPNN."""

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure


PEPTIDE_LENGTHS = (8, 10, 12, 14, 16, 18)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--experiment-dir",
        type=Path,
        required=True,
        help="Root directory of the monomer self-consistency experiment.",
    )
    parser.add_argument(
        "--run-name",
        default="run_01",
        help="LigandMPNN run directory to prepare (default: run_01).",
    )
    parser.add_argument(
        "--expected-count",
        type=int,
        default=200,
        help="Expected number of RFD3 backbones at each length.",
    )
    return parser.parse_args()


def model_number(path: Path) -> int:
    """Return the integer model number at the end of an RFD3 filename."""

    match = re.search(r"_model_(\d+)\.cif\.gz$", path.name)
    if match is None:
        raise ValueError(f"Could not read the model number from {path.name}")
    return int(match.group(1))


def read_cif(path: Path):
    """Read a gzip-compressed RFD3 CIF file with Biotite."""

    with gzip.open(path, mode="rt") as file:
        return get_structure(CIFFile.read(file), model=1)


def write_pdb(atom_array, path: Path) -> None:
    """Write an atom array as a PDB file for LigandMPNN."""

    pdb_file = PDBFile()
    pdb_file.set_structure(atom_array)
    pdb_file.write(path)


def prepare_length(
    experiment_dir: Path,
    run_name: str,
    peptide_length: int,
    expected_count: int,
) -> None:
    length_name = f"length_{peptide_length:02d}"
    source_dir = experiment_dir / "rfd3_backbones" / length_name
    run_dir = experiment_dir / "ligandmpnn_sequences" / run_name / length_name
    pdb_dir = run_dir / "input_pdbs"
    path_list = run_dir / "pdb_paths.json"

    # Calls model_number()
    source_paths = sorted(source_dir.glob("*.cif.gz"), key=model_number)
    if len(source_paths) != expected_count:
        raise ValueError(
            f"Expected {expected_count} structures in {source_dir}, "
            f"but found {len(source_paths)}"
        )

    pdb_dir.mkdir(parents=True, exist_ok=True)
    ligandmpnn_inputs = {}

    for source_path in source_paths:
        pdb_name = source_path.name.removesuffix(".cif.gz") + ".pdb"
        pdb_path = pdb_dir / pdb_name

        if not pdb_path.exists():
            atom_array = read_cif(source_path)
            ca_count = np.count_nonzero(atom_array.atom_name == "CA")
            if ca_count != peptide_length:
                raise ValueError(
                    f"Expected {peptide_length} CA atoms in {source_path}, "
                    f"but found {ca_count}"
                )
            write_pdb(atom_array, pdb_path)

        ligandmpnn_inputs[str(pdb_path.resolve())] = ""

    with path_list.open("w") as file:
        json.dump(ligandmpnn_inputs, file, indent=4)
        file.write("\n")

    print(f"Prepared {len(ligandmpnn_inputs)} inputs for {length_name}")
    print(f"Path list: {path_list}")


def main() -> None:
    args = parse_arguments()
    experiment_dir = args.experiment_dir.resolve()

    for peptide_length in PEPTIDE_LENGTHS:
        prepare_length(
            experiment_dir=experiment_dir,
            run_name=args.run_name,
            peptide_length=peptide_length,
            expected_count=args.expected_count,
        )


if __name__ == "__main__":
    main()
