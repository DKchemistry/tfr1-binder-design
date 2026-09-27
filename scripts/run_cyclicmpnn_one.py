#!/usr/bin/env python3

"""Run CyclicMPNN, and optionally Rosetta relaxation, for one backbone."""

import argparse
import csv
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def run_command(command, quiet=False):
    command = [str(item) for item in command]
    print(f"\n$ {shlex.join(command)}")

    if not quiet:
        subprocess.run(command, check=True)
        return

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, command)

    for line in result.stdout.splitlines():
        if line.startswith("Wrote "):
            print(line)


def get_selected_site(scores_csv, pdb_name):
    with open(scores_csv, newline="") as handle:
        rows = list(csv.DictReader(handle))

    selected = [
        row
        for row in rows
        if row["pdb"] == pdb_name
        and row["selected"].lower() in {"true", "1", "yes"}
    ]

    if len(selected) != 1:
        raise ValueError(
            f"Expected exactly one selected site for {pdb_name}; "
            f"found {len(selected)}"
        )

    row = selected[0]
    return row["chain"], int(row["sequence_index"])


def get_chain_ids(pdb_path):
    chain_ids = []

    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            chain = line[21].strip()
            if chain and chain not in chain_ids:
                chain_ids.append(chain)

    return chain_ids


def write_omit_json(
    output_path,
    structure_name,
    pdb_path,
    peptide_chain,
    site,
    allowed_aas,
):
    omitted_aas = "".join(
        amino_acid
        for amino_acid in AMINO_ACIDS
        if amino_acid not in allowed_aas
    )
    constraints = {
        chain: []
        for chain in get_chain_ids(pdb_path)
    }
    constraints[peptide_chain] = [
        [[site], omitted_aas]
    ]

    with open(output_path, "w") as handle:
        json.dump({structure_name: constraints}, handle)
        handle.write("\n")


def read_designed_sequence(fasta_path):
    sequences = []

    with open(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                sequences.append(line)

    if len(sequences) < 2:
        raise ValueError(f"Could not find designed sequence in {fasta_path}")

    sequence = sequences[-1]
    if "/" in sequence:
        raise ValueError("Unexpected multichain sequence in CyclicMPNN output.")

    return sequence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--scores", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--xml", type=Path)
    parser.add_argument("--design-chain", default="A")
    parser.add_argument(
        "--cyclicmpnn-dir",
        type=Path,
        default=Path("~/CyclicMPNN"),
    )
    parser.add_argument(
        "--model-name",
        default="cyclicmpnn_48_010",
    )
    parser.add_argument(
        "--thread-script",
        type=Path,
        default=SCRIPT_DIR / "thread_sequence.py",
    )
    parser.add_argument(
        "--relax-script",
        type=Path,
        default=SCRIPT_DIR / "run_relax.py",
    )
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

    args.cyclicmpnn_dir = args.cyclicmpnn_dir.expanduser().resolve()

    if args.no_relax and args.rounds != 1:
        raise ValueError("--no-relax requires --rounds 1")
    if not args.no_relax and args.xml is None:
        raise ValueError("--xml is required unless --no-relax is used")

    weights_dir = args.cyclicmpnn_dir / "cyclicmpnn_weights"
    weights_file = weights_dir / f"{args.model_name}.pt"
    if not weights_file.is_file():
        raise FileNotFoundError(f"CyclicMPNN weights not found: {weights_file}")

    if args.scores is None:
        chain = args.design_chain
        site = None
    else:
        chain, site = get_selected_site(args.scores, args.pdb.name)

    print(f"Input: {args.pdb}")
    print(f"CyclicMPNN model: {weights_file}")
    if site is None:
        print("Distal-site constraint: disabled")
    else:
        print(f"Selected lariat site: {chain}{site}")
        print(f"Allowed residues: {args.allowed_aas}")
    print(f"Rounds: {args.rounds}")
    print(f"Rosetta relaxation: {not args.no_relax}")

    args.output.mkdir(parents=True, exist_ok=True)
    current_input = args.pdb
    summary = []

    for round_number in range(1, args.rounds + 1):
        print(f"\n{'=' * 60}")
        print(f"ROUND {round_number}")
        print(f"{'=' * 60}")

        round_dir = args.output / f"round_{round_number}"
        if round_dir.exists() and any(round_dir.iterdir()):
            if args.force:
                shutil.rmtree(round_dir)
            else:
                raise SystemExit(
                    f"{round_dir} already contains files. "
                    "Use --force to overwrite."
                )
        round_dir.mkdir(parents=True, exist_ok=True)

        omit_json = None
        if site is not None:
            omit_json = round_dir / "omit_AA.jsonl"
            write_omit_json(
                output_path=omit_json,
                structure_name=current_input.stem,
                pdb_path=current_input,
                peptide_chain=chain,
                site=site,
                allowed_aas=args.allowed_aas,
            )

        command = [
            "conda", "run", "--no-capture-output",
            "-n", args.cyclicmpnn_env,
            "python",
            str(args.cyclicmpnn_dir / "protein_mpnn_run.py"),
            "--pdb_path", str(current_input),
            "--pdb_path_chains", chain,
            "--out_folder", str(round_dir),
            "--path_to_model_weights", str(weights_dir),
            "--model_name", args.model_name,
            "--num_seq_per_target", "1",
            "--batch_size", "1",
            "--sampling_temp", str(args.temperature),
            "--seed", str(args.seed),
        ]
        if omit_json is not None:
            command.extend(["--omit_AA_jsonl", str(omit_json)])

        run_command(command)

        fasta_path = round_dir / "seqs" / f"{current_input.stem}.fa"
        sequence = read_designed_sequence(fasta_path)
        print(f"\nRound {round_number} sequence:")
        print(sequence)

        selected_aa = ""
        if site is not None:
            selected_aa = sequence[site - 1]
            if selected_aa not in args.allowed_aas:
                raise ValueError(
                    f"{chain}{site} is {selected_aa}, "
                    f"but allowed residues are {args.allowed_aas}"
                )
            print(f"Lariat site {chain}{site}: {selected_aa} ✓")

        threaded_pdb = ""
        relaxed_pdb = ""

        if not args.no_relax:
            threaded_pdb = round_dir / "threaded.pdb"
            run_command([
                "conda", "run", "--no-capture-output",
                "-n", args.pyrosetta_env,
                "python", str(args.thread_script),
                str(current_input), str(threaded_pdb),
                "--chain", chain,
                "--sequence", sequence,
            ])

            relaxed_pdb = round_dir / "relaxed.pdb"
            print("\nRunning Rosetta Relax...")
            run_command([
                "conda", "run",
                "-n", args.pyrosetta_env,
                "python", str(args.relax_script),
                str(threaded_pdb), str(relaxed_pdb),
                "--xml", str(args.xml),
            ], quiet=True)

        summary.append({
            "round": round_number,
            "input_pdb": str(current_input),
            "sequence": sequence,
            "lariat_site": f"{chain}{site}" if site is not None else "",
            "lariat_residue": selected_aa,
            "threaded_pdb": str(threaded_pdb),
            "relaxed_pdb": str(relaxed_pdb),
        })

        if not args.no_relax:
            current_input = relaxed_pdb

    summary_path = args.output / "summary.tsv"
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=summary[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary)

    print(f"\n{'=' * 60}")
    print("COMPLETE")
    print(f"{'=' * 60}")
    print(f"Final sequence: {summary[-1]['sequence']}")
    if not args.no_relax:
        print(f"Final structure: {current_input}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
