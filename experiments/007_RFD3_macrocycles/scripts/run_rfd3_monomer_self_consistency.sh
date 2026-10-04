#!/usr/bin/env bash

set -euo pipefail

# Generate 200 cyclic monomer backbones at each peptide length used by the
# RFpeptides self-consistency benchmark.

FOUNDRY_DIR=/home/dkouv/foundry
RFD3_PYTHON="$FOUNDRY_DIR/.venv/bin/python"
RFD3_SCRIPT="$FOUNDRY_DIR/models/rfd3/src/rfd3/run_inference.py"
RFD3_CHECKPOINT=/home/dkouv/.foundry/checkpoints/rfd3_latest.ckpt
EXPECTED_FOUNDRY_COMMIT=69c70a88e0ed7348602ecfcef55d9ee6082cd675

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)/monomer_self_consistency"
CONFIG_DIR="$EXPERIMENT_DIR/configs/rfd3"
OUTPUT_ROOT="$EXPERIMENT_DIR/rfd3_backbones"

PEPTIDE_LENGTHS=(08 10 12 14 16 18)
DESIGNS_PER_LENGTH=200

actual_commit="$(git -C "$FOUNDRY_DIR" rev-parse HEAD)"
if [[ "$actual_commit" != "$EXPECTED_FOUNDRY_COMMIT" ]]; then
    echo "Foundry is checked out at the wrong commit."
    echo "Expected: $EXPECTED_FOUNDRY_COMMIT"
    echo "Actual:   $actual_commit"
    exit 1
fi

for peptide_length in "${PEPTIDE_LENGTHS[@]}"; do
    input_json="$CONFIG_DIR/length_${peptide_length}.json"
    output_dir="$OUTPUT_ROOT/length_${peptide_length}"

    mkdir -p "$output_dir"

    if find "$output_dir" -maxdepth 1 -name '*.cif.gz' -print -quit | grep -q .; then
        echo "Output structures already exist in: $output_dir"
        echo "Stopping rather than overwriting or silently skipping them."
        exit 1
    fi

    echo
    echo "Generating $DESIGNS_PER_LENGTH designs of length $peptide_length"
    echo "Output: $output_dir"

    "$RFD3_PYTHON" "$RFD3_SCRIPT" \
        out_dir="$output_dir" \
        ckpt_path="$RFD3_CHECKPOINT" \
        inputs="$input_json" \
        diffusion_batch_size="$DESIGNS_PER_LENGTH" \
        n_batches=1
done

echo
echo "RFD3 backbone generation complete."
