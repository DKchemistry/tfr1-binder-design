#!/usr/bin/env bash

set -euo pipefail

# Reproduce the AfCycDesign prediction conditions used for cyclic monomers:
# all five AlphaFold2 pTM models, 15% random masking, no dropout, and the
# type-2 cyclic positional offset. ColabDesign interprets six recycles as
# six recycled passes after the first forward pass (seven total passes).

AFCYC_PYTHON=/home/dkouv/miniforge3/envs/afcyc/bin/python
ALPHAFOLD_PARAMS=/home/dkouv/alphafold

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)/monomer_self_consistency"
PREDICTION_SCRIPT="$SCRIPT_DIR/run_afcyc_monomer_predictions.py"

INPUT_ROOT="$EXPERIMENT_DIR/ligandmpnn_sequences/run_01"
OUTPUT_ROOT="$EXPERIMENT_DIR/afcyc_predictions/run_01"

RECYCLES=6
RANDOM_MASK_FRACTION=0.15
SEED=0

export JAX_COMPILATION_CACHE_DIR="$HOME/.cache/jax"

exec "$AFCYC_PYTHON" "$PREDICTION_SCRIPT" \
    --input-root "$INPUT_ROOT" \
    --output-root "$OUTPUT_ROOT" \
    --params "$ALPHAFOLD_PARAMS" \
    --recycles "$RECYCLES" \
    --random-mask-fraction "$RANDOM_MASK_FRACTION" \
    --seed "$SEED" \
    "$@"
