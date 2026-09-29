# Experiment 006: Scaling the beta-paired 14-mer workflow

This experiment repeats the tensor generation, RFdiffusion, distal-site, and
ProteinMPNN settings from experiment 004 for 100 backbones. It intentionally
disables sequence threading and RosettaRelax. AfCyc therefore evaluates the
first ProteinMPNN sequence assignment in `round_1`.

## Configuration

The pipeline configuration is:

```text
experiments/006_scale_bpair_14mer/config/scale_1.toml
```

ProteinMPNN only needs the repository root:

```toml
[paths]
proteinmpnn = "/home/dkouv/ProteinMPNN"
```

The pipeline calls `protein_mpnn_run.py` below that directory. ProteinMPNN then
loads its default `v_48_020.pt` checkpoint from its own
`vanilla_model_weights/` directory. A separate weights path is not needed for
this experiment.

The no-relax settings are:

```toml
[mpnn_relax]
method = "proteinmpnn"
rounds = 1
relax = false
```

`rounds` must be `1` when relaxation is disabled because later rounds in the
iterative workflow normally consume the previous round's relaxed structure.

## Run command

From the repository root:

```sh
python3 scripts/run_pipeline.py \
  --config experiments/006_scale_bpair_14mer/config/scale_1.toml
```

The run will write stage-specific outputs beneath:

```text
experiments/006_scale_bpair_14mer/outputs/scale-1/
├── tensors/
├── rfdiffusion/
├── distal_site/
├── mpnn_relax/
│   └── <design>/
│       ├── round_1/
│       │   ├── omit_AA.jsonl
│       │   └── seqs/<design>.fa
│       └── summary.tsv
└── oracle/
    └── <design>/round_1/
```

There should be no `threaded.pdb` or `relaxed.pdb` files in this experiment.
Generated files below `outputs/` are ignored by Git.
