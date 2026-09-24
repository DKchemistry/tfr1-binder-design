#!/usr/bin/env python3

import argparse
import csv
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def run_command(command, quiet=False):
    print(f"\n$ {shlex.join(str(x) for x in command)}")

    if not quiet:
        subprocess.run(command, check=True)
        return

    # Rosetta is extremely verbose. Hide successful output, but show it on error.
    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise subprocess.CalledProcessError(
            result.returncode,
            command,
        )

    # Keep our useful final line from run_relax.py.
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
        aa for aa in AMINO_ACIDS
        if aa not in allowed_aas
    )

    constraints = {
        chain: []
        for chain in get_chain_ids(pdb_path)
    }

    constraints[peptide_chain] = [
        [[site], omitted_aas]
    ]

    data = {
        structure_name: constraints
    }

    with open(output_path, "w") as handle:
        json.dump(data, handle)
        handle.write("\n")


def read_mpnn_sequence(fasta_path):
    sequences = []

    with open(fasta_path) as handle:
        for line in handle:
            line = line.strip()

            if line and not line.startswith(">"):
                sequences.append(line)

    if len(sequences) < 2:
        raise ValueError(
            f"Could not find designed sequence in {fasta_path}"
        )

    # ProteinMPNN writes the input sequence first,
    # followed by generated sequence(s).
    sequence = sequences[-1]

    if "/" in sequence:
        raise ValueError(
            "Unexpected multichain sequence in ProteinMPNN output."
        )

    return sequence


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--xml", type=Path, required=True)

    parser.add_argument(
        "--proteinmpnn-dir",
        type=Path,
        default=Path("~/work/ProteinMPNN"),
    )

    parser.add_argument(
        "--thread-script",
        type=Path,
        default=Path("scripts/thread_sequence.py"),
    )

    parser.add_argument(
        "--relax-script",
        type=Path,
        default=Path("scripts/run_relax.py"),
    )

    parser.add_argument("--rounds", type=int, default=4)
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=1)

    parser.add_argument(
        "--allowed-aas",
        default="CDEK",
        help="Allowed residues at lariat attachment site",
    )

    parser.add_argument(
        "--mpnn-env",
        default="proteinmpnn",
    )

    parser.add_argument(
        "--pyrosetta-env",
        default="pyrosetta",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing round directories",
    )

    args = parser.parse_args()

    args.proteinmpnn_dir = args.proteinmpnn_dir.expanduser()

    chain, site = get_selected_site(
        args.scores,
        args.pdb.name,
    )

    print(f"Input: {args.pdb}")
    print(f"Selected lariat site: {chain}{site}")
    print(f"Allowed residues: {args.allowed_aas}")
    print(f"Rounds: {args.rounds}")

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

        # ------------------------------------------------------------
        # 1. Write the residue-specific ProteinMPNN constraint
        # ------------------------------------------------------------

        omit_json = round_dir / "omit_AA.jsonl"

        write_omit_json(
            output_path=omit_json,
            structure_name=current_input.stem,
            pdb_path=current_input,
            peptide_chain=chain,
            site=site,
            allowed_aas=args.allowed_aas,
        )

        # ------------------------------------------------------------
        # 2. ProteinMPNN
        # ------------------------------------------------------------

        mpnn_command = [
            "conda", "run", "--no-capture-output",
            "-n", args.mpnn_env,
            "python",
            str(args.proteinmpnn_dir / "protein_mpnn_run.py"),

            "--pdb_path", str(current_input),
            "--pdb_path_chains", chain,

            "--omit_AA_jsonl", str(omit_json),

            "--out_folder", str(round_dir),

            "--num_seq_per_target", "1",
            "--batch_size", "1",

            "--sampling_temp", str(args.temperature),
            "--seed", str(args.seed),
        ]

        run_command(mpnn_command)

        fasta_path = (
            round_dir
            / "seqs"
            / f"{current_input.stem}.fa"
        )

        sequence = read_mpnn_sequence(fasta_path)

        print(f"\nRound {round_number} sequence:")
        print(sequence)

        selected_aa = sequence[site - 1]

        if selected_aa not in args.allowed_aas:
            raise ValueError(
                f"{chain}{site} is {selected_aa}, "
                f"but allowed residues are {args.allowed_aas}"
            )

        print(
            f"Lariat site {chain}{site}: "
            f"{selected_aa} ✓"
        )

        # ------------------------------------------------------------
        # 3. Thread the sequence
        # ------------------------------------------------------------

        threaded_pdb = round_dir / "threaded.pdb"

        thread_command = [
            "conda", "run", "--no-capture-output",
            "-n", args.pyrosetta_env,
            "python",
            str(args.thread_script),

            str(current_input),
            str(threaded_pdb),

            "--chain", chain,
            "--sequence", sequence,
        ]

        run_command(thread_command)

        # ------------------------------------------------------------
        # 4. Rosetta FastRelax
        # ------------------------------------------------------------

        relaxed_pdb = round_dir / "relaxed.pdb"

        relax_command = [
            "conda", "run",
            "-n", args.pyrosetta_env,
            "python",
            str(args.relax_script),

            str(threaded_pdb),
            str(relaxed_pdb),

            "--xml", str(args.xml),
        ]

        print("\nRunning Rosetta Relax...")
        run_command(relax_command, quiet=True)

        summary.append(
            {
                "round": round_number,
                "input_pdb": str(current_input),
                "sequence": sequence,
                "lariat_site": f"{chain}{site}",
                "lariat_residue": selected_aa,
                "threaded_pdb": str(threaded_pdb),
                "relaxed_pdb": str(relaxed_pdb),
            }
        )

        # Relaxed structure becomes next round's input.
        current_input = relaxed_pdb

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

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

    for row in summary:
        print(
            f"Round {row['round']}: "
            f"{row['sequence']} "
            f"({row['lariat_site']}={row['lariat_residue']})"
        )

    print(f"\nFinal structure: {current_input}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()