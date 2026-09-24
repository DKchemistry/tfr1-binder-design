#!/usr/bin/env python3

import argparse

import pyrosetta
from pyrosetta.toolbox.mutants import mutate_residue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdb")
    parser.add_argument("output_pdb")
    parser.add_argument("--chain", required=True)
    parser.add_argument("--sequence", required=True)
    args = parser.parse_args()

    pyrosetta.init("-mute all")

    pose = pyrosetta.pose_from_pdb(args.input_pdb)

    # Find all protein residues belonging to the requested PDB chain.
    chain_positions = [
        i
        for i in range(1, pose.total_residue() + 1)
        if pose.pdb_info().chain(i) == args.chain
        and pose.residue(i).is_protein()
    ]

    if len(chain_positions) != len(args.sequence):
        raise ValueError(
            f"Chain {args.chain} has {len(chain_positions)} residues, "
            f"but sequence has {len(args.sequence)} residues."
        )

    print(f"Threading chain {args.chain}")
    print(f"Sequence: {args.sequence}")

    for pose_index, amino_acid in zip(chain_positions, args.sequence):
        pdb_number = pose.pdb_info().number(pose_index)

        old_aa = pose.residue(pose_index).name1()

        print(
            f"{args.chain}{pdb_number}: "
            f"{old_aa} -> {amino_acid}"
        )

        if old_aa != amino_acid:
            mutate_residue(
                pose,
                pose_index,
                amino_acid,
                pack_radius=0.0,
            )

    pose.dump_pdb(args.output_pdb)

    print(f"Wrote {args.output_pdb}")


if __name__ == "__main__":
    main()