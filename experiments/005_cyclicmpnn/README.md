I set up cyclicmpnn like this:

# Clone your known-working ProteinMPNN environment
conda create -n cyclicmpnn --clone proteinmpnn -y

# Activate the separate CyclicMPNN environment
conda activate cyclicmpnn

# Clone CyclicMPNN directly into your home directory
cd ~
git clone https://github.com/ParisaH-Lab/CyclicMPNN.git
cd ~/CyclicMPNN

# Verify GPU/PyTorch
python -c "import torch; print('PyTorch:', torch.__version__); print('CUDA:', torch.version.cuda); print('CUDA available:', torch.cuda.is_available()); print('GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'None')"

# Verify the script runs
python protein_mpnn_run.py --help

Everything looks good, lets check the weights:

cd ~/CyclicMPNN
ls -lh cyclicmpnn_weights/

total 20M
-rw-r--r-- 1 dkouv dkouv 20M Sep 26 16:50 cyclicmpnn_48_010.pt

Looks good too.

Now, we will write the first exp to: `experiments/005_cyclicmpnn/outputs/cycmpnn-1`.

We need to carefully modify our `scripts/run_pipeline.py` to be able to run analyses like this that we can specify in a TOML that remain backwards compatabile. The script must still be easy to read. We have diffusion outputs already here: `experiments/004_oracle_success/outputs/exp1/rfdiffusion`, we need to chain that to predict a backbone with cyclicmpnn this time, we will not be running RosettaRelax, we will simply use the 1 backbone assignment and feed it directly to the oracle. We want to perserve as much similar directory structure as we can, so we can still make a round_1 in the what will be the mpnn_relax folder, we will just not proceed with rosetta relax.

## Implemented workflow

The pipeline configuration for this run is:

```text
experiments/005_cyclicmpnn/configs/cycmpnn_1_pipeline.toml
```

Run it from the repository root with:

```sh
python scripts/run_pipeline.py \
  --config experiments/005_cyclicmpnn/configs/cycmpnn_1_pipeline.toml
```

This configuration skips tensor generation, RFdiffusion, and distal-site
selection. It reads the experiment-004 RFdiffusion PDBs and distal-site CSV,
runs one CyclicMPNN sequence assignment per backbone using
`cyclicmpnn_48_010.pt`, skips threading and Rosetta relaxation, and sends each
`round_1` sequence to AfCyc.

The new results are written below:

```text
experiments/005_cyclicmpnn/outputs/cycmpnn-1/
├── mpnn_relax/
│   └── <design>/
│       ├── round_1/
│       │   ├── omit_AA.jsonl
│       │   └── seqs/<design>.fa
│       └── summary.tsv
└── oracle/
    └── <design>/round_1/
```

There is deliberately no `threaded.pdb` or `relaxed.pdb` in this run. AfCyc
reads the designed sequence from `summary.tsv`, so neither file is required.

## Compatibility notes

Existing pipeline TOMLs remain valid. When `[stages]` is absent, all stages
run. When `method` and `relax` are absent from `[mpnn_relax]`, they default to
`proteinmpnn` and `true`, respectively. Consequently the experiment-004 TOML
continues to request its original four ProteinMPNN/RosettaRelax rounds.

When relaxation is disabled, `rounds` must be `1`: later rounds in the original
workflow depend on the previous round's relaxed backbone.

## Run result: cycmpnn-1

The configured pipeline completed successfully on 2026-09-27.

- CyclicMPNN sequence assignments: 100/100
- AfCyc predictions: 100/100
- Sequence length: 14 residues for every design
- Distal-site constraint: satisfied for every design
- Threaded structures: 0, as configured
- Relaxed structures: 0, as configured
- AfCyc artifacts: 100 each of `prediction.pdb`, `metrics.json`, `pae.npy`,
  `cyclic_offset.npy`, and `.complete`

The sequence-design stage finished at approximately 05:56 CEST and the AfCyc
stage ran from approximately 05:56 to 06:38 CEST. As a quick output-integrity
summary, normalized iPAE ranged down to 0.4250 (best: `rfdiffusion_46`) with a
median of 0.8785; binder pLDDT had a median of 0.6217 and maximum of 0.9191.
These unfiltered values only confirm readable oracle metrics. Structural RMSD
and interface analysis should be performed separately before judging designs.


## Back to me 

Seems iPAE did not fair as well as before. 

```sh
conda run -n biotite python scripts/analyze_rmsd_ipae.py \
  --oracle experiments/005_cyclicmpnn/outputs/cycmpnn-1/oracle \
  --reference experiments/004_oracle_success/outputs/exp1/rfdiffusion \
  --oracle-binder-ch B \
  --oracle-target-ch A \
  --reference-binder-ch A \
  --reference-target-ch B \
  --round 1 \
  --dpi 300 \
  --hotspot-residues 20-25 \
  --output-dir experiments/005_cyclicmpnn/outputs/cycmpnn-1/rmsd_ipae_analysis
```

