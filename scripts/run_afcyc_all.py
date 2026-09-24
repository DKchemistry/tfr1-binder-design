#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def get_round_sequence(summary_tsv, round_number=None):
    with open(summary_tsv, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        raise ValueError(
            f"No rounds found in {summary_tsv}"
        )

    if round_number is None:
        row = max(
            rows,
            key=lambda row: int(row["round"]),
        )
    else:
        matches = [
            row for row in rows
            if int(row["round"]) == round_number
        ]

        if len(matches) != 1:
            raise ValueError(
                f"Expected one round-{round_number} row in "
                f"{summary_tsv}; found {len(matches)}"
            )

        row = matches[0]

    return int(row["round"]), row["sequence"]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--design-dir",
        type=Path,
        required=True,
        help="Directory containing iterative-design subdirectories",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--params",
        type=Path,
        default=Path("~/alphafold"),
    )

    parser.add_argument(
        "--afcyc-script",
        type=Path,
        default=SCRIPT_DIR / "afcyc_predict.py",
    )

    parser.add_argument(
        "--round",
        type=int,
        default=None,
        help="MPNN round to evaluate. Default: latest completed round.",
    )

    parser.add_argument(
        "--target-pdb",
        type=Path,
        required=True,
        help="PDB containing the target/receptor structure for AfCyc.",
    )

    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--recycles", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--afcyc-env", default="afcyc")

    args = parser.parse_args()

    args.target_pdb = args.target_pdb.expanduser().resolve()

    if not args.target_pdb.is_file():
        raise FileNotFoundError(
            f"Target PDB not found: {args.target_pdb}"
        )

    args.params = args.params.expanduser()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    design_dirs = sorted(
        path
        for path in args.design_dir.iterdir()
        if path.is_dir()
    )

    env = os.environ.copy()
    env["JAX_COMPILATION_CACHE_DIR"] = str(
        Path("~/.cache/jax").expanduser()
    )

    failures = []

    for design_dir in design_dirs:
        name = design_dir.name

        summary_tsv = design_dir / "summary.tsv"

        if not summary_tsv.exists():
            print(f"Skipping {name}: missing summary")
            continue

        round_number, sequence = get_round_sequence(
            summary_tsv,
            args.round,
        )

        output = (
            args.output_dir
            / name
            / f"round_{round_number}"
        )
        complete = output / ".complete"
        log_file = output / "run.log"

        if complete.exists():
            print(f"Skipping {name}: already completed")
            continue

        output.mkdir(parents=True, exist_ok=True)

        print("\n" + "=" * 70)
        print(f"AfCyc: {name}")
        print(f"Sequence: {sequence}")
        print(f"Length: {len(sequence)}")
        print(f"Target structure: {args.target_pdb}")
        print("=" * 70)

        command = [
            "conda", "run", "--no-capture-output",
            "-n", args.afcyc_env,
            "python",
            str(args.afcyc_script),

            "--pdb", str(args.target_pdb),
            "--target_chain", args.target_chain,
            "--sequence", sequence,
            "--params", str(args.params),
            "--out", str(output),
            "--recycles", str(args.recycles),
            "--seed", str(args.seed),
        ]

        with open(log_file, "w") as log:
            result = subprocess.run(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
            )

        if result.returncode == 0:
            complete.touch()
            print(f"Completed: {name}")
        else:
            failures.append(name)
            print(f"FAILED: {name}")
            print(f"See: {log_file}")

    print("\n" + "=" * 70)

    if failures:
        print("Finished with failures:")
        for name in failures:
            print(f"  {name}")
    else:
        print("All AfCyc predictions completed successfully.")


if __name__ == "__main__":
    main()