#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from math import pi
from math import radians
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import rotaxis
from scipy.spatial import cKDTree


def heavy_atoms(residue):
    """Yield non-hydrogen atoms."""
    for atom in residue.get_atoms():
        element = (atom.element or "").strip().upper()
        if element != "H" and not atom.get_name().upper().startswith("H"):
            yield atom


def unit(vector):
    """Return a unit vector."""
    return vector / np.linalg.norm(vector)


def sidechain_direction(residue):
    """
    Return the unit vector pointing from CA toward the side-chain side.

    If a real CB exists, use CA -> CB.

    For glycine, construct a virtual CB direction by rotating CA -> N
    by -120 degrees around the CA -> C axis. This is the same
    pseudo-CB construction used by Bio.PDB for glycine exposure.
    """
    ca = residue["CA"].get_vector()

    if "CB" in residue:
        vector = (residue["CB"].get_vector() - ca).get_array()
        return unit(vector)

    n_from_ca = residue["N"].get_vector() - ca
    c_from_ca = residue["C"].get_vector() - ca

    rotation = rotaxis(radians(-120.0), c_from_ca)
    virtual_cb_vector = n_from_ca.left_multiply(rotation)

    return unit(virtual_cb_vector.get_array())


def peptide_self_clearance(probe, peptide_residues, residue_index):
    """
    Distance from the probe to non-local peptide atoms.

    Ignore the residue itself and its immediate cyclic neighbours.
    """
    n_residues = len(peptide_residues)

    excluded = {
        residue_index,
        (residue_index - 1) % n_residues,
        (residue_index + 1) % n_residues,
    }

    coords = []

    for i, residue in enumerate(peptide_residues):
        if i in excluded:
            continue

        coords.extend(atom.coord for atom in heavy_atoms(residue))

    if not coords:
        return float("inf")

    coords = np.asarray(coords)

    return float(
        np.min(np.linalg.norm(coords - probe, axis=1))
    )


def score_structure(pdb_path, peptide_chain_id="A", probe_distance=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, pdb_path)
    model = next(structure.get_models())

    peptide_chain = model[peptide_chain_id]

    peptide_residues = []

    for residue in peptide_chain:

        is_standard_residue = residue.id[0] == " "

        has_backbone_atoms = (
            "N" in residue
            and "CA" in residue
            and "C" in residue
        )

        if is_standard_residue and has_backbone_atoms:
            peptide_residues.append(residue)
            
    # Everything outside the peptide chain is treated as receptor.
    receptor_atoms = []

    for chain in model:
        if chain.id == peptide_chain_id:
            continue

        for residue in chain:
            receptor_atoms.extend(heavy_atoms(residue))

    receptor_coords = np.asarray(
        [atom.coord for atom in receptor_atoms]
    )

    receptor_tree = cKDTree(receptor_coords)

    rows = []

    for index, residue in enumerate(peptide_residues):
        direction = sidechain_direction(residue)

        ca = np.asarray(residue["CA"].coord)

        # Generic point representing the direction in which
        # a future side chain / derivatization handle would extend.
        probe = ca + probe_distance * direction

        ca_clearance, nearest_index = receptor_tree.query(ca)
        probe_clearance, _ = receptor_tree.query(probe)

        # Vector from CA toward the nearest receptor atom.
        receptor_direction = unit(
            receptor_coords[nearest_index] - ca
        )

        # +1 = directly away from receptor
        #  0 = tangential
        # -1 = directly toward receptor
        outward_cosine = -float(
            np.dot(direction, receptor_direction)
        )

        self_clearance = peptide_self_clearance(
            probe,
            peptide_residues,
            index,
        )

        rows.append(
            {
                "pdb": pdb_path.name,
                "chain": peptide_chain_id,
                "residue": residue.id[1],
                "sequence_index": index + 1,
                "receptor_CA_clearance_A": float(ca_clearance),
                "receptor_probe_clearance_A": float(probe_clearance),
                "clearance_gain_A": float(
                    probe_clearance - ca_clearance
                ),
                "outward_cosine": outward_cosine,
                "peptide_probe_clearance_A": self_clearance,
                "probe_x": float(probe[0]),
                "probe_y": float(probe[1]),
                "probe_z": float(probe[2]),
            }
        )

    return rows


def choose_site(rows, clearance_tolerance=0.5):
    """
    Select one site.

    First maximize receptor clearance.

    Positions within clearance_tolerance Angstrom of the best
    clearance are treated as effectively tied. Among those,
    prefer:

      1. most outward-pointing side-chain direction
      2. greatest clearance from the rest of the peptide
      3. lowest sequence index
    """
    best_clearance = max(
        row["receptor_probe_clearance_A"]
        for row in rows
    )

    finalists = [
        row
        for row in rows
        if (
            best_clearance
            - row["receptor_probe_clearance_A"]
            <= clearance_tolerance
        )
    ]

    return max(
        finalists,
        key=lambda row: (
            row["outward_cosine"],
            row["peptide_probe_clearance_A"],
            -row["sequence_index"],
        ),
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "pdb_dir",
        type=Path,
        help="Directory containing PDB files",
    )

    parser.add_argument(
        "--peptide-chain",
        default="A",
        help="Peptide chain ID (default: A)",
    )

    parser.add_argument(
        "--probe-distance",
        type=float,
        default=4.0,
        help="Generic side-chain probe distance in Angstrom (default: 4.0)",
    )

    parser.add_argument(
        "--clearance-tolerance",
        type=float,
        default=0.5,
        help="Clearance difference considered a near-tie (default: 0.5 A)",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional CSV containing scores for every residue",
    )

    args = parser.parse_args()

    pdb_paths = sorted(args.pdb_dir.glob("*.pdb"))

    if not pdb_paths:
        raise SystemExit(
            f"No PDB files found in {args.pdb_dir}"
        )

    all_rows = []

    for pdb_path in pdb_paths:
        rows = score_structure(
            pdb_path,
            peptide_chain_id=args.peptide_chain,
            probe_distance=args.probe_distance,
        )

        selected = choose_site(
            rows,
            clearance_tolerance=args.clearance_tolerance,
        )

        for row in rows:
            row["selected"] = row is selected
            all_rows.append(row)

        print(
            f"{pdb_path.name}: "
            f"{args.peptide_chain}{selected['residue']} "
            f"(sequence index {selected['sequence_index']}, "
            f"probe clearance {selected['receptor_probe_clearance_A']:.2f} Å, "
            f"outward cosine {selected['outward_cosine']:.2f}, "
            f"peptide clearance {selected['peptide_probe_clearance_A']:.2f} Å, "
            f"probe xyz = "
            f"({selected['probe_x']:.3f}, "
            f"{selected['probe_y']:.3f}, "
            f"{selected['probe_z']:.3f}))"
        )
        
    if args.output:
        args.output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=all_rows[0].keys(),
            )
            writer.writeheader()
            writer.writerows(all_rows)


if __name__ == "__main__":
    main()