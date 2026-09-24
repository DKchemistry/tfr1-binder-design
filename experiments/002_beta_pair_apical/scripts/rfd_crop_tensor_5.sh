#!/usr/bin/env bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv

mkdir -p "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_5"

cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_5/rfd_crop_tensor_5" \
  inference.num_designs=100 \
  inference.input_pdb="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  'contigmap.contigs=[14-14 A189-383/0]' \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A210-213_bl-14_bsl-4" \
  scaffoldguided.scaffold_list="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A210-213_bl-14_bsl-4.txt" \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A209,A210,A211,A212]'