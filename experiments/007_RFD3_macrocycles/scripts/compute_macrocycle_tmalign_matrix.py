#!/usr/bin/env python3

"""Compute the all-by-all TM-align matrix for the 1,200 RFD3 backbones.

For each pair of structures, the matrix entry is the average of the two
TM-scores reported by TM-align: one normalized by each structure's length.
The resulting matrix is symmetric and has ones on its diagonal.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import subprocess
from pathlib import Path

import biotite.structure as struc
import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure


PEPTIDE_LENGTHS = (8, 10, 12, 14, 16, 18)
EXPECTED_BACKBONES_PER_LENGTH = 200

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1] / "monomer_self_consistency"
DEFAULT_BACKBONE_DIR = EXPERIMENT_DIR / "rfd3_backbones"
DEFAULT_OUTPUT_DIR = EXPERIMENT_DIR / "structural_tsne"
DEFAULT_TMALIGN = Path(
    "/home/dkouv/miniforge3/envs/macrocycle-tsne/bin/TMalign"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--backbone-dir",
        type=Path,
        default=DEFAULT_BACKBONE_DIR,
        help=f"RFD3 backbone directory. Default: {DEFAULT_BACKBONE_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Structural t-SNE directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--tmalign",
        type=Path,
        default=DEFAULT_TMALIGN,
        help=f"TM-align executable. Default: {DEFAULT_TMALIGN}",
    )
    return parser.parse_args()


def model_number(path: Path) -> int:
    """Return the integer model number at the end of an RFD3 filename."""

    match = re.search(r"_model_(\d+)\.cif\.gz$", path.name)
    if match is None:
        raise ValueError(f"Could not read model number from {path.name}")
    return int(match.group(1))


def find_backbones(backbone_dir: Path) -> list[tuple[int, Path]]:
    """Find all backbones in peptide-length and model-number order."""

    backbones = []

    for peptide_length in PEPTIDE_LENGTHS:
        length_dir = backbone_dir / f"length_{peptide_length:02d}"
        paths = sorted(length_dir.glob("*.cif.gz"), key=model_number)

        if len(paths) != EXPECTED_BACKBONES_PER_LENGTH:
            raise ValueError(
                f"Expected {EXPECTED_BACKBONES_PER_LENGTH} backbones in "
                f"{length_dir}; found {len(paths)}"
            )

        backbones.extend((peptide_length, path) for path in paths)

    return backbones


def read_cif(path: Path) -> struc.AtomArray:
    """Read the first model from a gzip-compressed CIF file."""

    with gzip.open(path, mode="rt") as handle:
        return get_structure(CIFFile.read(handle), model=1)


def write_ca_pdb(source_path: Path, output_path: Path, peptide_length: int) -> None:
    """Write the C-alpha atoms used by TM-align to a small PDB file."""

    structure = read_cif(source_path)
    ca_atoms = structure[structure.atom_name == "CA"]

    if len(ca_atoms) != peptide_length:
        raise ValueError(
            f"Expected {peptide_length} C-alpha atoms in {source_path}; "
            f"found {len(ca_atoms)}"
        )

    pdb_file = PDBFile()
    pdb_file.set_structure(ca_atoms)
    pdb_file.write(output_path)


def prepare_tmalign_inputs(
    backbones: list[tuple[int, Path]],
    output_dir: Path,
) -> tuple[list[dict[str, object]], Path, Path]:
    """Create standardized PDB inputs, an index, and a TM-align file list."""

    input_dir = output_dir / "tmalign_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    index_rows = []

    for matrix_index, (peptide_length, source_path) in enumerate(backbones):
        backbone = source_path.name.removesuffix(".cif.gz")
        pdb_path = input_dir / f"{backbone}.pdb"
        write_ca_pdb(source_path, pdb_path, peptide_length)

        index_rows.append(
            {
                "matrix_index": matrix_index,
                "peptide_length": peptide_length,
                "backbone": backbone,
                "source_path": str(source_path.resolve()),
                "tmalign_pdb_path": str(pdb_path.resolve()),
            }
        )

    index_path = output_dir / "structure_index.tsv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(index_rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(index_rows)

    list_path = output_dir / "tmalign_input_files.txt"
    with list_path.open("w") as handle:
        for row in index_rows:
            handle.write(f"{Path(str(row['tmalign_pdb_path'])).name}\n")

    return index_rows, input_dir, list_path


def run_tmalign(
    executable: Path,
    input_dir: Path,
    list_path: Path,
    raw_output_path: Path,
) -> None:
    """Run TM-align's all-against-all mode and save its compact output."""

    command = [
        str(executable),
        "-dir",
        f"{input_dir.resolve()}/",
        str(list_path.resolve()),
        "-outfmt",
        "2",
    ]

    print("Running TM-align all-against-all comparison...")
    with raw_output_path.open("w") as output_handle:
        subprocess.run(
            command,
            stdout=output_handle,
            check=True,
            text=True,
        )


def parse_tmalign_output(
    raw_output_path: Path,
    index_rows: list[dict[str, object]],
    output_dir: Path,
) -> np.ndarray:
    """Average the two reported TM-scores and construct a symmetric matrix."""

    name_to_index = {
        Path(str(row["tmalign_pdb_path"])).name: int(row["matrix_index"])
        for row in index_rows
    }
    index_to_backbone = {
        int(row["matrix_index"]): str(row["backbone"])
        for row in index_rows
    }

    structure_count = len(index_rows)
    matrix = np.full((structure_count, structure_count), np.nan, dtype=float)
    np.fill_diagonal(matrix, 1.0)

    pair_rows = []
    seen_pairs = set()

    with raw_output_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) != 11:
                raise ValueError(f"Unexpected TM-align output line:\n{line}")

            name_i = Path(fields[0]).name
            name_j = Path(fields[1]).name
            index_i = name_to_index[name_i]
            index_j = name_to_index[name_j]

            if index_i == index_j:
                continue

            pair = tuple(sorted((index_i, index_j)))
            tm_score_i = float(fields[2])
            tm_score_j = float(fields[3])
            mean_tm_score = (tm_score_i + tm_score_j) / 2.0

            if pair in seen_pairs:
                if not np.isclose(matrix[index_i, index_j], mean_tm_score):
                    raise ValueError(f"Conflicting scores for matrix pair {pair}")
                continue

            seen_pairs.add(pair)
            matrix[index_i, index_j] = mean_tm_score
            matrix[index_j, index_i] = mean_tm_score

            pair_rows.append(
                {
                    "matrix_index_i": index_i,
                    "matrix_index_j": index_j,
                    "backbone_i": index_to_backbone[index_i],
                    "backbone_j": index_to_backbone[index_j],
                    "tm_score_i": tm_score_i,
                    "tm_score_j": tm_score_j,
                    "mean_tm_score": mean_tm_score,
                }
            )

    expected_pairs = structure_count * (structure_count - 1) // 2
    if len(seen_pairs) != expected_pairs:
        raise ValueError(
            f"Expected {expected_pairs:,} unique structure pairs; "
            f"found {len(seen_pairs):,}"
        )
    if np.isnan(matrix).any():
        raise ValueError("The completed TM-score matrix contains missing values")
    if not np.allclose(matrix, matrix.T):
        raise ValueError("The completed TM-score matrix is not symmetric")
    if not np.allclose(np.diag(matrix), 1.0):
        raise ValueError("The completed TM-score matrix diagonal is not one")
    if matrix.min() < 0.0 or matrix.max() > 1.0:
        raise ValueError("The completed TM-score matrix contains values outside [0, 1]")

    pair_path = output_dir / "tmalign" / "pair_scores.tsv"
    pair_path.parent.mkdir(parents=True, exist_ok=True)
    with pair_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(pair_rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(pair_rows)

    return matrix


def write_metadata(
    path: Path,
    executable: Path,
    structure_count: int,
    matrix: np.ndarray,
) -> None:
    """Record the matrix definition and basic validation results."""

    version_result = subprocess.run(
        [str(executable), "-v"],
        check=False,
        capture_output=True,
        text=True,
    )
    version = (version_result.stdout or version_result.stderr).strip()

    metadata = {
        "structure_count": structure_count,
        "unique_pair_count": structure_count * (structure_count - 1) // 2,
        "matrix_definition": (
            "mean of TM-align scores normalized by each structure length"
        ),
        "symmetric": bool(np.allclose(matrix, matrix.T)),
        "diagonal": 1.0,
        "minimum_tm_score": float(matrix.min()),
        "maximum_tm_score": float(matrix.max()),
        "tmalign_executable": str(executable.resolve()),
        "tmalign_version": version,
    }

    with path.open("w") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")


def main() -> None:
    args = parse_arguments()
    backbone_dir = args.backbone_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    executable = args.tmalign.expanduser().resolve()

    if not executable.is_file():
        raise FileNotFoundError(f"Could not find TM-align executable:\n{executable}")

    output_dir.mkdir(parents=True, exist_ok=True)
    backbones = find_backbones(backbone_dir)
    index_rows, input_dir, list_path = prepare_tmalign_inputs(
        backbones,
        output_dir,
    )

    tmalign_dir = output_dir / "tmalign"
    tmalign_dir.mkdir(parents=True, exist_ok=True)
    raw_output_path = tmalign_dir / "raw_tmalign_output.tsv"
    run_tmalign(executable, input_dir, list_path, raw_output_path)

    matrix = parse_tmalign_output(raw_output_path, index_rows, output_dir)
    matrix_path = tmalign_dir / "tm_scores.npy"
    np.save(matrix_path, matrix)

    metadata_path = tmalign_dir / "metadata.json"
    write_metadata(metadata_path, executable, len(backbones), matrix)

    print(f"Structures: {len(backbones):,}")
    print(f"Unique pairs: {len(backbones) * (len(backbones) - 1) // 2:,}")
    print(f"TM-score range: {matrix.min():.4f} to {matrix.max():.4f}")
    print(f"Saved: {matrix_path}")
    print(f"Saved: {metadata_path}")


if __name__ == "__main__":
    main()
