#!/usr/bin/env python3

import argparse
import json
import re
import shutil
from importlib.metadata import version
from pathlib import Path


PEPTIDE_LENGTHS = (8, 10, 12, 14, 16, 18)
MODEL_NAMES = tuple(f"model_{number}_ptm" for number in range(1, 6))
SEQUENCE_ID_PATTERN = re.compile(r"(?:^|,\s*)id=(\d+)(?:,|$)")


def read_fasta(fasta_path):
    """Return the header and sequence for every record in a FASTA file."""

    records = []
    header = None
    sequence_lines = []

    with fasta_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_lines)))

                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)

    if header is not None:
        records.append((header, "".join(sequence_lines)))

    return records


def read_designed_sequences(fasta_path, expected_length):
    """Read the eight LigandMPNN designs, excluding the native record."""

    designs = {}

    for header, sequence in read_fasta(fasta_path):
        match = SEQUENCE_ID_PATTERN.search(header)

        # The first FASTA record is the input sequence and has no id field.
        if match is None:
            continue

        sequence_id = int(match.group(1))
        designs[sequence_id] = sequence.upper()

    expected_ids = set(range(1, 9))

    if set(designs) != expected_ids:
        raise ValueError(
            f"Expected sequence IDs 1-8 in {fasta_path}; "
            f"found {sorted(designs)}"
        )

    wrong_lengths = {
        sequence_id: len(sequence)
        for sequence_id, sequence in designs.items()
        if len(sequence) != expected_length
    }

    if wrong_lengths:
        observed = ", ".join(
            f"id={sequence_id}: {sequence_length}"
            for sequence_id, sequence_length in sorted(wrong_lengths.items())
        )
        raise ValueError(
            f"Expected length {expected_length} in {fasta_path}; "
            f"found {observed}"
        )

    return designs


def collect_prediction_tasks(input_root, output_root):
    """Collect valid sequences in a stable length/file/sequence order."""

    tasks = []
    excluded_backbones = []

    for peptide_length in PEPTIDE_LENGTHS:
        length_name = f"length_{peptide_length:02d}"
        fasta_dir = input_root / length_name / "seqs"

        if not fasta_dir.is_dir():
            raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")

        fasta_paths = sorted(fasta_dir.glob("*.fa"))

        if len(fasta_paths) != 200:
            raise ValueError(
                f"Expected 200 FASTA files in {fasta_dir}; "
                f"found {len(fasta_paths)}"
            )

        for fasta_path in fasta_paths:
            try:
                designs = read_designed_sequences(
                    fasta_path,
                    expected_length=peptide_length,
                )
            except ValueError as error:
                excluded_backbones.append((fasta_path.stem, str(error)))
                continue

            for sequence_id, sequence in sorted(designs.items()):
                output_dir = (
                    output_root
                    / length_name
                    / fasta_path.stem
                    / f"sequence_{sequence_id:02d}"
                )

                tasks.append(
                    {
                        "length": peptide_length,
                        "backbone": fasta_path.stem,
                        "fasta_path": fasta_path,
                        "sequence_id": sequence_id,
                        "sequence": sequence,
                        "output_dir": output_dir,
                    }
                )

    return tasks, excluded_backbones


def add_cyclic_offset(model):
    """Apply the type-2 cyclic relative positional encoding."""

    import numpy as np

    peptide_length = model._len
    residue_numbers = np.arange(peptide_length)

    repeated_numbers = np.stack(
        [residue_numbers, residue_numbers + peptide_length],
        axis=-1,
    )

    linear_offset = (
        residue_numbers[:, None]
        - residue_numbers[None, :]
    )

    cyclic_distance = np.abs(
        repeated_numbers[:, None, :, None]
        - repeated_numbers[None, :, None, :]
    ).min(axis=(2, 3))

    wrapped_is_shorter = cyclic_distance < np.abs(linear_offset)
    cyclic_distance[wrapped_is_shorter] *= -1

    model._inputs["offset"] = (
        cyclic_distance * np.sign(linear_offset)
    )


def install_jax_compatibility_shims():
    """Support the newer JAX version in the existing AfCyc environment."""

    import jax
    import jax.numpy as jnp

    if not hasattr(jax, "tree_map"):
        jax.tree_map = jax.tree_util.tree_map

    original_clip = jnp.clip

    def clip_compat(
        arr=None,
        a_min=None,
        a_max=None,
        *,
        min=None,
        max=None,
    ):
        lower = min if min is not None else a_min
        upper = max if max is not None else a_max
        return original_clip(arr, min=lower, max=upper)

    jnp.clip = clip_compat

    return jax


def predict_sequence(model, task, args):
    """Run all five AlphaFold2 pTM models for one peptide sequence."""

    import numpy as np

    output_dir = task["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Reset to the same seed for each sequence. Successive model calls then
    # receive deterministic, but distinct, random keys for random masking.
    model.set_seed(args.seed)

    model_metrics = []

    for model_name in MODEL_NAMES:
        aux = model.predict(
            seq=task["sequence"],
            models=[model_name],
            num_recycles=args.recycles,
            sample_models=False,
            dropout=False,
            return_aux=True,
            verbose=False,
        )

        prediction_path = output_dir / f"{model_name}.pdb"
        pae_path = output_dir / f"{model_name}_pae.npy"

        model.save_current_pdb(str(prediction_path))

        pae = np.asarray(aux["pae"])
        np.save(pae_path, pae)

        model_metrics.append(
            {
                "model": model_name,
                "plddt": float(aux["log"]["plddt"]),
                "plddt_100": float(aux["log"]["plddt"] * 100.0),
                "mean_pae_angstrom": float(pae.mean()),
                "ptm": float(aux["log"]["ptm"]),
            }
        )

    best_model = max(model_metrics, key=lambda item: item["plddt"])
    best_model_name = best_model["model"]

    shutil.copyfile(
        output_dir / f"{best_model_name}.pdb",
        output_dir / "prediction.pdb",
    )
    shutil.copyfile(
        output_dir / f"{best_model_name}_pae.npy",
        output_dir / "pae.npy",
    )

    np.save(
        output_dir / "cyclic_offset.npy",
        np.asarray(model._inputs["offset"]),
    )

    metrics = {
        "backbone": task["backbone"],
        "source_fasta": str(task["fasta_path"].resolve()),
        "sequence_id": task["sequence_id"],
        "sequence": task["sequence"],
        "peptide_length": task["length"],
        "cyclic_offset_type": 2,
        "models": list(MODEL_NAMES),
        "best_model_by_plddt": best_model_name,
        "best_plddt": best_model["plddt"],
        "best_plddt_100": best_model["plddt_100"],
        "recycles": args.recycles,
        "total_network_passes": args.recycles + 1,
        "random_mask_fraction": args.random_mask_fraction,
        "dropout": False,
        "seed": args.seed,
        "colabdesign_version": version("colabdesign"),
        "model_metrics": model_metrics,
    }

    with (output_dir / "metrics.json").open("w") as handle:
        json.dump(metrics, handle, indent=2)

    # This marker is written last. Its presence means the prediction is safe
    # to skip when the campaign is restarted.
    (output_dir / ".complete").touch()

    return best_model


def print_task_summary(tasks, excluded_backbones):
    completed = sum(
        (task["output_dir"] / ".complete").is_file()
        for task in tasks
    )

    print(f"Valid predictions:     {len(tasks)}")
    print(f"Already complete:      {completed}")
    print(f"Remaining predictions: {len(tasks) - completed}")
    print(f"Excluded backbones:    {len(excluded_backbones)}")

    for backbone, reason in excluded_backbones:
        print(f"  {backbone}")
        print(f"    {reason}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Predict LigandMPNN-designed cyclic monomers with AfCycDesign."
        )
    )

    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--params", type=Path, required=True)
    parser.add_argument("--recycles", type=int, default=6)
    parser.add_argument("--random-mask-fraction", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--max-predictions",
        type=int,
        default=None,
        help="Stop after this many new predictions. Useful for a test run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and report progress without loading AfCyc.",
    )

    args = parser.parse_args()

    args.input_root = args.input_root.expanduser().resolve()
    args.output_root = args.output_root.expanduser().resolve()
    args.params = args.params.expanduser().resolve()

    if args.recycles < 0:
        raise ValueError("--recycles must be at least zero")

    if not 0.0 <= args.random_mask_fraction <= 1.0:
        raise ValueError("--random-mask-fraction must be between zero and one")

    tasks, excluded_backbones = collect_prediction_tasks(
        args.input_root,
        args.output_root,
    )

    print_task_summary(tasks, excluded_backbones)

    if args.dry_run:
        return

    if not args.params.is_dir():
        raise FileNotFoundError(f"AlphaFold data directory not found: {args.params}")

    if not (args.params / "params").is_dir():
        raise FileNotFoundError(
            f"AlphaFold parameter directory not found: {args.params / 'params'}"
        )

    jax = install_jax_compatibility_shims()

    gpu_devices = [
        device
        for device in jax.devices()
        if device.platform == "gpu"
    ]

    if not gpu_devices:
        raise RuntimeError("AfCyc cannot see a GPU. No predictions were started.")

    print(f"GPU: {gpu_devices[0]}", flush=True)
    print(
        f"Recycles: {args.recycles} "
        f"({args.recycles + 1} total network passes)",
        flush=True,
    )
    print(f"Models: {', '.join(MODEL_NAMES)}", flush=True)
    print(f"Random masking: {args.random_mask_fraction:.0%}", flush=True)

    # Import ColabDesign only after the inexpensive input validation above.
    from colabdesign import mk_afdesign_model

    model = mk_afdesign_model(
        protocol="hallucination",
        data_dir=str(args.params),
        model_names=list(MODEL_NAMES),
        use_mlm=True,
    )
    model.set_opt(mlm_dropout=args.random_mask_fraction)

    current_length = None
    new_predictions = 0

    try:
        for task_number, task in enumerate(tasks, start=1):
            complete_marker = task["output_dir"] / ".complete"

            if complete_marker.is_file():
                continue

            if (
                args.max_predictions is not None
                and new_predictions >= args.max_predictions
            ):
                break

            if task["length"] != current_length:
                current_length = task["length"]
                model.prep_inputs(length=current_length)
                add_cyclic_offset(model)

            print(
                f"[{task_number}/{len(tasks)}] "
                f"{task['backbone']} sequence {task['sequence_id']}",
                flush=True,
            )

            best_model = predict_sequence(model, task, args)
            new_predictions += 1

            print(
                f"  complete: {best_model['model']}, "
                f"pLDDT={best_model['plddt']:.4f}",
                flush=True,
            )

    except KeyboardInterrupt:
        print(
            "\nStopped by user. Restart the same command to resume.",
            flush=True,
        )
        return

    print(f"New predictions completed: {new_predictions}", flush=True)
    print_task_summary(tasks, excluded_backbones)


if __name__ == "__main__":
    main()
