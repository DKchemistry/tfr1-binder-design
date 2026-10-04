#!/usr/bin/env bash

set -euo pipefail

# Reproduce the LigandMPNN conditions used for the RFpeptides monomer
# self-consistency benchmark: eight sequences per backbone at temperature 0.1.

LIGANDMPNN_DIR=/home/dkouv/LigandMPNN
LIGANDMPNN_PYTHON="$LIGANDMPNN_DIR/.venv/bin/python"
LIGANDMPNN_SCRIPT="$LIGANDMPNN_DIR/run.py"
LIGANDMPNN_CHECKPOINT="$LIGANDMPNN_DIR/model_params/ligandmpnn_v_32_010_25.pt"
EXPECTED_LIGANDMPNN_COMMIT=26ec57ac976ade5379920dbd43c7f97a91cf82de

BIOTITE_PYTHON=/home/dkouv/miniforge3/envs/biotite/bin/python

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)/monomer_self_consistency"
PREPARE_SCRIPT="$SCRIPT_DIR/prepare_ligandmpnn_inputs.py"
OUTPUT_ROOT="$EXPERIMENT_DIR/ligandmpnn_sequences/run_01"

PEPTIDE_LENGTHS=(08 10 12 14 16 18)
BACKBONES_PER_LENGTH=200
SEQUENCES_PER_BACKBONE=8
TEMPERATURE=0.1
SEED=2026

actual_commit="$(git -C "$LIGANDMPNN_DIR" rev-parse HEAD)"
if [[ "$actual_commit" != "$EXPECTED_LIGANDMPNN_COMMIT" ]]; then
    echo "LigandMPNN is checked out at the wrong commit."
    echo "Expected: $EXPECTED_LIGANDMPNN_COMMIT"
    echo "Actual:   $actual_commit"
    exit 1
fi

cuda_device="$($LIGANDMPNN_PYTHON -c 'import torch; print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")')"
if [[ -z "$cuda_device" ]]; then
    echo "LigandMPNN cannot see a CUDA device. Stopping before sequence design."
    exit 1
fi
echo "CUDA device: $cuda_device"

"$BIOTITE_PYTHON" "$PREPARE_SCRIPT" \
    --experiment-dir "$EXPERIMENT_DIR" \
    --run-name run_01 \
    --expected-count "$BACKBONES_PER_LENGTH"

for peptide_length in "${PEPTIDE_LENGTHS[@]}"; do
    output_dir="$OUTPUT_ROOT/length_${peptide_length}"
    pdb_path_list="$output_dir/pdb_paths.json"

    if compgen -G "$output_dir/seqs/*.fa" > /dev/null; then
        echo "LigandMPNN sequence outputs already exist in: $output_dir/seqs"
        echo "Stopping rather than overwriting them."
        exit 1
    fi

    echo
    echo "Designing $SEQUENCES_PER_BACKBONE sequences for each length-$peptide_length backbone"
    echo "Output: $output_dir"

    "$LIGANDMPNN_PYTHON" "$LIGANDMPNN_SCRIPT" \
        --model_type ligand_mpnn \
        --checkpoint_ligand_mpnn "$LIGANDMPNN_CHECKPOINT" \
        --pdb_path_multi "$pdb_path_list" \
        --out_folder "$output_dir" \
        --temperature "$TEMPERATURE" \
        --batch_size "$SEQUENCES_PER_BACKBONE" \
        --number_of_batches 1 \
        --seed "$SEED"

    fasta_count=$(find "$output_dir/seqs" -maxdepth 1 -name '*.fa' | wc -l)
    if [[ "$fasta_count" -ne "$BACKBONES_PER_LENGTH" ]]; then
        echo "Expected $BACKBONES_PER_LENGTH FASTA files, but found $fasta_count."
        exit 1
    fi
done

echo
echo "LigandMPNN sequence design complete."
