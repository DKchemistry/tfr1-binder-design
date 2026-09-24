#!/usr/bin/env bash

set -euo pipefail

ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion
PROTEINMPNN=/home/dkouv/ProteinMPNN
AF_PARAMS=/home/dkouv/alphafold

PROJECT="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/test-simple-pipeline"

TARGET_PDB="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb"

TENSOR_DIR="$PROJECT/tensors/interface_tensors_dir"
TENSOR_LIST="$PROJECT/tensors/interface_tensors_dir.txt"

mkdir -p \
  "$PROJECT/tensors" \
  "$PROJECT/rfdiffusion" \
  "$PROJECT/distal_site" \
  "$PROJECT/mpnn_relax" \
  "$PROJECT/oracle"


echo
echo "======================================================================"
echo "1. Interface tensors"
echo "======================================================================"

conda run --no-capture-output -n pyrosetta \
  python "$ROOT/scripts/interface_tensors/make_interface_tensor.py" \
  --input_pdb "$TARGET_PDB" \
  --out_dir "$TENSOR_DIR" \
  --binderlen 14 \
  --target_adj A210-213 \
  --binder_ss E \
  --binder_ss_len 4


echo
echo "======================================================================"
echo "2. RFdiffusion"
echo "======================================================================"

(
  cd "$RFD"

  conda run --no-capture-output -n SE3nv \
    python scripts/run_inference.py \
    --config-name base \
    inference.output_prefix="$PROJECT/rfdiffusion/rfdiffusion" \
    inference.num_designs=5 \
    inference.input_pdb="$TARGET_PDB" \
    'contigmap.contigs=[14-14 A189-383/0]' \
    inference.model_runner=ScaffoldedSampler \
    scaffoldguided.scaffoldguided=True \
    scaffoldguided.scaffold_dir="$TENSOR_DIR" \
    scaffoldguided.scaffold_list="$TENSOR_LIST" \
    inference.cyclic=True \
    inference.cyc_chains=a \
    diffuser.T=50 \
    'ppi.hotspot_res=[A209,A210,A211,A212]'
)


echo
echo "======================================================================"
echo "3. Distal-site selection"
echo "======================================================================"

conda run --no-capture-output -n biopython \
  python "$ROOT/scripts/select_distal_site.py" \
  "$PROJECT/rfdiffusion" \
  --peptide-chain A \
  --output "$PROJECT/distal_site/distal_site.csv"


echo
echo "======================================================================"
echo "4. ProteinMPNN / RosettaRelax"
echo "======================================================================"

conda run --no-capture-output -n biopython \
  python "$ROOT/scripts/run_mpnn_relax_all.py" \
  --pdb-dir "$PROJECT/rfdiffusion" \
  --scores "$PROJECT/distal_site/distal_site.csv" \
  --output-dir "$PROJECT/mpnn_relax" \
  --xml "$ROOT/RosettaRelax/fast_relax_binder_A.xml" \
  --proteinmpnn-dir "$PROTEINMPNN" \
  --rounds 4


echo
echo "======================================================================"
echo "5. AfCyc oracle"
echo "======================================================================"

XLA_PYTHON_CLIENT_PREALLOCATE=false \
conda run --no-capture-output -n biopython \
  python "$ROOT/scripts/run_afcyc_all.py" \
  --design-dir "$PROJECT/mpnn_relax" \
  --output-dir "$PROJECT/oracle" \
  --target-pdb "$TARGET_PDB" \
  --params "$AF_PARAMS" \
  --afcyc-env afcyc \
  --recycles 5


echo
echo "======================================================================"
echo "Pipeline complete"
echo "======================================================================"
echo "Project: $PROJECT"