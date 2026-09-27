#!/usr/bin/env python3

"""Run the single-backbone CyclicMPNN workflow over a PDB directory."""

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir", type=Path, required=True)
    parser.add_argument("--scores", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--xml", type=Path)
    parser.add_argument("--cyclicmpnn-dir", type=Path, required=True)
    parser.add_argument("--design-chain", default="A")
    parser.add_argument("--model-name", default="cyclicmpnn_48_010")
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--allowed-aas", default="CDEK")
    parser.add_argument("--cyclicmpnn-env", default="cyclicmpnn")
    parser.add_argument("--pyrosetta-env", default="pyrosetta")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--no-relax",
        action="store_true",
        help="Generate one sequence without threading or Rosetta relaxation.",
    )
    args = parser.parse_args()

    if args.no_relax and args.rounds != 1:
        raise ValueError("--no-relax requires --rounds 1")
    if not args.no_relax and args.xml is None:
        raise ValueError("--xml is required unless --no-relax is used")

    pdb_paths = sorted(args.pdb_dir.glob("*.pdb"))
    if not pdb_paths:
        raise SystemExit(f"No PDB files found in {args.pdb_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for pdb_path in pdb_paths:
        design_name = pdb_path.stem
        output_path = args.output_dir / design_name
        summary_path = output_path / "summary.tsv"
        final_pdb = output_path / f"round_{args.rounds}" / "relaxed.pdb"
        is_complete = (
            summary_path.exists()
            if args.no_relax
            else final_pdb.exists()
        )

        if is_complete and not args.force:
            print(f"Skipping {design_name}: requested output already exists")
            continue

        print(f"\n{'#' * 70}")
        print(f"Running {design_name}")
        print(f"{'#' * 70}\n")

        command = [
            sys.executable,
            str(SCRIPT_DIR / "run_cyclicmpnn_one.py"),
            "--pdb", str(pdb_path),
            "--output", str(output_path),
            "--cyclicmpnn-dir", str(args.cyclicmpnn_dir),
            "--design-chain", args.design_chain,
            "--model-name", args.model_name,
            "--rounds", str(args.rounds),
            "--temperature", str(args.temperature),
            "--seed", str(args.seed),
            "--allowed-aas", args.allowed_aas,
            "--cyclicmpnn-env", args.cyclicmpnn_env,
            "--pyrosetta-env", args.pyrosetta_env,
        ]
        if args.scores is not None:
            command.extend(["--scores", str(args.scores)])
        if args.xml is not None:
            command.extend(["--xml", str(args.xml)])
        if args.no_relax:
            command.append("--no-relax")
        if args.force:
            command.append("--force")

        subprocess.run(command, check=True)

    print("\nAll available designs completed.")


if __name__ == "__main__":
    main()
