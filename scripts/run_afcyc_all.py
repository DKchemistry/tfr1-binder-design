#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
from pathlib import Path


def get_round4_sequence(summary_tsv):
    with open(summary_tsv, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    round4 = [
        row for row in rows
        if int(row["round"]) == 4
    ]

    if len(round4) != 1:
        raise ValueError(
            f"Expected one round-4 row in {summary_tsv}; "
            f"found {len(round4)}"
        )

    return round4[0]["sequence"]


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
        "--afcyc-script",
        type=Path,
        default=Path("scripts/afcyc_predict.py"),
    )

    parser.add_argument(
        "--params",
        type=Path,
        default=Path("~/alphafold"),
    )

    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--recycles", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--afcyc-env", default="afcyc")

    args = parser.parse_args()

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
        pdb = design_dir / "round_4" / "relaxed.pdb"

        if not summary_tsv.exists() or not pdb.exists():
            print(f"Skipping {name}: missing summary or final PDB")
            continue

        sequence = get_round4_sequence(summary_tsv)

        output = args.output_dir / name
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
        print(f"Structure: {pdb}")
        print("=" * 70)

        command = [
            "conda", "run", "--no-capture-output",
            "-n", args.afcyc_env,
            "python",
            str(args.afcyc_script),

            "--pdb", str(pdb),
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