#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--pdb-dir", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--xml", type=Path, required=True)

    args = parser.parse_args()

    pdb_paths = sorted(args.pdb_dir.glob("*.pdb"))

    if not pdb_paths:
        raise SystemExit(
            f"No PDB files found in {args.pdb_dir}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for pdb_path in pdb_paths:
        design_name = pdb_path.stem
        output_path = args.output_dir / design_name

        final_pdb = output_path / "round_4" / "relaxed.pdb"

        if final_pdb.exists():
            print(
                f"Skipping {design_name}: "
                "round 4 relaxed structure already exists"
            )
            continue

        print(f"\n{'#' * 70}")
        print(f"Running {design_name}")
        print(f"{'#' * 70}\n")

        command = [
            "python",
            "scripts/run_mpnn_relax_one.py",
            "--pdb",
            str(pdb_path),
            "--scores",
            str(args.scores),
            "--output",
            str(output_path),
            "--xml",
            str(args.xml),
        ]

        subprocess.run(
            command,
            check=True,
        )

    print("\nAll available designs completed.")


if __name__ == "__main__":
    main()