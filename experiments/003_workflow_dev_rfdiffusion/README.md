I need to make a re-useable workflow that goes through the following process: 

# 1. SS/ADJ Tensor Calculations for Target X and Binder Y.

This should be optional in the event I want to forego beta pairing. 

Otherwise, the related arguements are the same as before. I'll start running this workflow again here (in `003_workflow_dev_rfdiffusion`)

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate pyrosetta

python "$ROOT/scripts/interface_tensors/make_interface_tensor.py" \
  --input_pdb "$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  --out_dir "$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/interface_tensors_A210-213_bl-14_bsl-4" \
  --binderlen 14 \
  --target_adj A210-213 \
  --binder_ss E \
  --binder_ss_len 4
```

We may be able to keep most of these argument names, except for `out_dir`, which should instead be `tensor_out_dir`. It also need to nested one level deeper than I expect. 

When we run: 

```sh
  --out_dir "$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/interface_tensors_A210-213_bl-14_bsl-4" \
```
We save one file a level above the `interface_tensors_A210-213_bl-14_bsl-4` dir. Like so, 

```
ls -l outputs/
interface_tensors_A210-213_bl-14_bsl-4/
interface_tensors_A210-213_bl-14_bsl-4.txt
```

So, `tensor_out_dir` should be a level lower than I expect. I want a run to be self-contained under one directory structure. Let's imagine the pipeline takes an argument for the project directory; we will call it `project_1` here.

If I want to get the tensors out, I need:

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate pyrosetta

python "$ROOT/scripts/interface_tensors/make_interface_tensor.py" \
  --input_pdb "$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  --out_dir "$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/project_1/tensors/interface_tensors_A210-213_bl-14_bsl-4" \
  --binderlen 14 \
  --target_adj A210-213 \
  --binder_ss E \
  --binder_ss_len 4
```

This gives:

```text
outputs/
└── project_1/
    └── tensors/
        ├── interface_tensors_A210-213_bl-14_bsl-4/
        │   ├── 6wrw_ChA_189-383_adj0_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj0_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj1_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj1_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj2_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj2_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj3_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj3_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj4_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj4_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj5_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj5_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj6_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj6_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj7_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj7_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj8_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj8_1_ss.pt
        │   ├── 6wrw_ChA_189-383_adj9_1_adj.pt
        │   └── 6wrw_ChA_189-383_adj9_1_ss.pt
        └── interface_tensors_A210-213_bl-14_bsl-4.txt
```

The important layout is therefore:

```text
outputs/
└── project_1/
    └── tensors/
        ├── <tensor_out_dir>/
        │   └── *.pt
        └── <tensor_out_dir>.txt
```

The `.txt` file is written alongside the tensor output directory, so `tensor_out_dir` itself needs to sit one level below `tensors/`.

# 2. RFdiffusion for Macrocycles / Scaffold Conditioning

Previously, I used a shell script like this to generate macrocycles of a fixed length, optionally using SS/ADJ conditioning tensors and specifying hotspot residues:

```sh
#!/usr/bin/env bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv

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
```

In the eventual pipeline, many of these RFdiffusion arguments should be configurable.

The SS/ADJ scaffold conditioning should be optional:

```sh
inference.model_runner=ScaffoldedSampler
scaffoldguided.scaffoldguided=True
scaffoldguided.scaffold_dir=<tensor_dir>
scaffoldguided.scaffold_list=<tensor_list>
```

However, `inference.input_pdb` should **not** be independently supplied to this stage. The pipeline should derive it from the project/run configuration so that tensor generation and RFdiffusion cannot accidentally operate on different input PDBs.

Likewise, RFdiffusion output paths should be opinionated rather than supplied directly by the user. A pipeline run should remain self-contained within the project directory rather than placing RFdiffusion outputs at the top level of that project directory.

For example:

```text
outputs/
└── project_1/
    ├── tensors/
    │   ├── interface_tensors_A210-213_bl-14_bsl-4/
    │   │   └── *.pt
    │   └── interface_tensors_A210-213_bl-14_bsl-4.txt
    │
    └── rfdiffusion/
        ├── <prefix>_0.pdb
        ├── <prefix>_0.trb
        ├── <prefix>_1.pdb
        ├── <prefix>_1.trb
        ├── ...
        └── traj/
```

Conceptually:

```text
project_dir/
├── tensors/
│   └── ...
└── rfdiffusion/
    └── RFdiffusion outputs
```

`project_dir` represents one run of the eventual pipeline. Since we only intend to try one RFdiffusion parameter set per pipeline run, there is no need for another RFdiffusion run directory beneath `rfdiffusion/`.

The pipeline should construct an `inference.output_prefix` equivalent to:

```text
<project_dir>/rfdiffusion/<prefix>
```

For example:

```text
inference.output_prefix=
    <project_dir>/rfdiffusion/rfdiffusion
```

which would produce files like:

```text
project_1/
└── rfdiffusion/
    ├── rfdiffusion_0.pdb
    ├── rfdiffusion_0.trb
    ├── rfdiffusion_1.pdb
    ├── rfdiffusion_1.trb
    ├── ...
    └── traj/
```

When scaffold conditioning is enabled, the corresponding tensor paths should similarly be derived from the same project directory:

```text
scaffoldguided.scaffold_dir=
    <project_dir>/tensors/<tensor_out_dir>

scaffoldguided.scaffold_list=
    <project_dir>/tensors/<tensor_out_dir>.txt
```

This keeps the relationship between the input PDB, generated tensors, and RFdiffusion outputs under control of the pipeline rather than requiring the caller to manually coordinate paths between stages.

The important distinction is that **`project_1` is the pipeline run**, while `tensors/` and `rfdiffusion/` are stage-specific outputs within that run.

# 3. Distal Site Selection 

In a previous project, I wrote a `select_distal_site.py` script to select a distal site that ProteinMPNN will require to be a C, D, E, or K. It operates on the PDB files generated by RFdiffusion. 

Previously, I ran it like so: 

```sh
conda activate biopython
python scripts/select_distal_site.py \
    proof-of-concept/rfd_tf1r_macrocycle/mps \
    --peptide-chain B \
    --output proof-of-concept/rfd_tf1r_macrocycle/mps/distal_site_scores.csv
```

`--peptide-chain` is the binder. Now, I will try it on some of the backbones I generated during `experiment_002`. 

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
conda activate biopython
python "$ROOT/scripts/select_distal_site.py" \
    "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_5" \
    --peptide-chain A \
    --output "$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/project_1/distal_site/distal_site_scores.csv"
```
It works great. This should also be opinionated. The pipeline now knows what the RFdiffusion output directory is. It also knows where the equivalent of `project_1` is, so it should just append `distal_site/distal_site_scores.csv`. 

The CSV itself looks like this: 

```text
pdb,chain,residue,sequence_index,receptor_CA_clearance_A,receptor_probe_clearance_A,clearance_gain_A,outward_cosine,peptide_probe_clearance_A,probe_x,probe_y,probe_z,selected
rfd_crop_tensor_5_0.pdb,A,1,1,10.96742893822457,14.789719306413375,3.8222903681888045,0.9397291278311708,6.913935898368953,24.53263495088201,-28.966479780177316,-0.930999426550098,True
rfd_crop_tensor_5_0.pdb,A,2,2,9.831269397561954,11.289710178626875,1.4584407810649207,0.1882221012783239,5.482898755867348,22.409848871823787,-27.337788683828997,6.868293184025108,False
```

Every residue of the PDB is evaluated and we only need to select the true condition to assign it. 

Running this should be optional, as not every run will require it.

# 3. ProteinMPNN and RosettaRelax

This part seems most complex.

If desired, the distal site is used to constrain that position to having a particular amino acid at that site using the `--omit_AA_jsonl` argument in proteinmpnn. 

If desired, the target and binder is threaded and RosettaRelax is used for relax the structure. This process may repeat in iterative fashion. 

In my previous project, a combination of scripts accomplished this and used an XML file to set up the Rosetta parameters. I have downloaded them here: 

```sh
RosettaRelax/fast_relax_binder_A.xml
scripts/thread_sequence.py
scripts/run_relax.py
scripts/run_mpnn_relax_one.py
scripts/run_mpnn_relax_all.py
```
These scripts only support some of proteinmpnn's potential arguments. It should be possible to supply the others, though I do not have the need for them now. Another design consideration is the temperature setting. At very low temperature sampling, proteinmpnn has an underflow condition. Using `--sampling_temp 0.1` has been in safe in my hands so far, but some papers report using very low temperatures that I have tried to incorporate. To prevent a design from failing, we would want to (1) allow dynamically retrying the sequence assignment with higher temperatures until it succeeds/an upper bound is reached, and/or (2) skipping backbones that can not have a sequence assigned to it at a given temperature. 

The other issue is handling the directory layout so everything is still under whatever the project dir was set as. 

The last proteinmpnn assignment should be carried into the oracle stage. 

# 4. AlphaFold2 Cyclic Oracle 

Once the final ProteinMPNN sequence has been assigned for the binder, we begin the cyclic oracle. 

This was previously handled by: 

```sh
scripts/run_afcyc_all.py
scripts/afcyc_predict.py
```

This is more straightforward because these settings are somewhat fixed. Though, there is the issue of what sequence we are going to try to fold. It is advised (typically) to fold the full chains/complex. I don't have the VRAM to do that with TfR1. So, we should be able to supply the sequence we want to use for prediction. 

# An example piecemeal workflow on this machine. 

To better understand the how the workflow above will work, let's start simple. 

First, I will just run the components piecemeal on a small scale and iterate on it until I get the desired directory structure out. 

## Paths I will need now or in the future: 

### Pre-installed environments. 

The workflow will need to use various environments, so we need to be able to set them in a config file. 

Environments: 

pyrosetta_env:`conda activate pyrosetta`
rfdifussion_env: `conda activate SE3nv`

### Scripts

Certain scripts will be called, so we need to be able to access them and so they should be in the config.

Scripts we cannot edit the internals of: 

interface_tensor_script: `/home/dkouv/work/tfr1-binder-design/scripts/interface_tensors/make_interface_tensor.py`

rfd_inference_script: `/home/dkouv/RFdiffusion/scripts/run_inference.py`

### Paths

A project directory dir, here it is: 

project_dir_path: `/home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project`

Pre-cleaned input pdb: `/home/dkouv/work/tfr1-binder-design/data/PDBs/6wrw_ChA_189-383.pdb`

## SS/ADJ Calculation 

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate pyrosetta

python "$ROOT/scripts/interface_tensors/make_interface_tensor.py" \
  --input_pdb "$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  --out_dir "$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/tensors/interface_tensors_dir" \
  --binderlen 14 \
  --target_adj A210-213 \
  --binder_ss E \
  --binder_ss_len 4
```
Gives: 

```
outputs/
└── piecemeal-project/
    └── tensors/
        ├── interface_tensors_dir/
        │   ├── *_adj.pt
        │   └── *_ss.pt
        └── interface_tensors_dir.txt
```
More concretely, the tensor stage produces:

```
<pipeline_run>/
└── tensors/
    ├── <tensor_out_dir>/
    │   ├── <design>_adj.pt
    │   ├── <design>_ss.pt
    │   └── ...
    └── <tensor_out_dir>.txt
```

For this run:

outputs/
└── piecemeal-project/
    └── tensors/
        ├── interface_tensors_dir/
        │   ├── 6wrw_ChA_189-383_adj0_1_adj.pt
        │   ├── 6wrw_ChA_189-383_adj0_1_ss.pt
        │   ├── ...
        │   ├── 6wrw_ChA_189-383_adj9_1_adj.pt
        │   └── 6wrw_ChA_189-383_adj9_1_ss.pt
        └── interface_tensors_dir.txt

The important point is that tensors/ contains both:

the tensor directory itself, containing the generated .pt files
a matching .txt file alongside that directory

## RFdiffusion 

We will start with 1 design to test pipeline execution.

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv

cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/rfdiffusion/rfdiffusion" \
  inference.num_designs=1 \
  inference.input_pdb="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  'contigmap.contigs=[14-14 A189-383/0]' \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/tensors/interface_tensors_dir" \
  scaffoldguided.scaffold_list="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/tensors/interface_tensors_dir.txt" \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A209,A210,A211,A212]'
```
## RFdiffusion

We will start with one design to test pipeline execution.

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv

cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/rfdiffusion/rfdiffusion" \
  inference.num_designs=1 \
  inference.input_pdb="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  'contigmap.contigs=[14-14 A189-383/0]' \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/tensors/interface_tensors_dir" \
  scaffoldguided.scaffold_list="$ROOT/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/tensors/interface_tensors_dir.txt" \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A209,A210,A211,A212]'
```

The RFdiffusion stage adds a new `rfdiffusion/` directory alongside the existing `tensors/` directory:

```text
piecemeal-project/
├── tensors/
│   └── ...                     # already covered
│
└── rfdiffusion/
    ├── rfdiffusion_0.pdb
    ├── rfdiffusion_0.trb
    └── traj/
        ├── rfdiffusion_0_Xt-1_traj.pdb
        └── rfdiffusion_0_pX0_traj.pdb
```

Conceptually, this stage produces:

```text
<pipeline_run>/
├── tensors/
│   └── ...
└── rfdiffusion/
    ├── <prefix>_<design>.pdb
    ├── <prefix>_<design>.trb
    └── traj/
        ├── <prefix>_<design>_Xt-1_traj.pdb
        └── <prefix>_<design>_pX0_traj.pdb
```

For this pipeline, `rfdiffusion/` is the RFdiffusion output directory for the run. There is no additional run-specific directory beneath it because the pipeline itself represents a single RFdiffusion parameter set.

The filename prefix is supplied through:

```text
inference.output_prefix=
    <pipeline_run>/rfdiffusion/<prefix>
```

For this test:

```text
inference.output_prefix=
    piecemeal-project/rfdiffusion/rfdiffusion
```

With one design, this produces `rfdiffusion_0.pdb` and `rfdiffusion_0.trb`, along with the corresponding trajectory files under `traj/`.

The existing tensor outputs are used as inputs to this stage but are unchanged, so their contents do not need to be repeated here.

## Distal Site 

I have a script that uses a simple heuristic to identify a residue site that will later be used as a synthetic handle.

The script assumes that its input directory contains a set of `*.pdb` files. In the eventual pipeline, this input directory will always be the RFdiffusion output directory.

```sh
conda activate biopython

cd /home/dkouv/work/tfr1-binder-design/

python scripts/select_distal_site.py \
    /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/rfdiffusion \
    --peptide-chain A \
    --output /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/distal_site/distal_site.csv
```

This stage reads the RFdiffusion-generated PDB files and adds a new `distal_site/` directory containing the selected residue information:

```text
piecemeal-project/
├── tensors/
│   └── ...                     # already covered
│
├── rfdiffusion/
│   └── ...                     # already covered
│
└── distal_site/
    └── distal_site.csv
```

Conceptually:

```text
<pipeline_run>/
├── tensors/
│   └── ...
├── rfdiffusion/
│   └── RFdiffusion outputs
└── distal_site/
    └── distal_site.csv
```

The input directory should be derived from the pipeline run directory:

```text
<pipeline_run>/rfdiffusion/
```

and the output path should likewise be opinionated:

```text
<pipeline_run>/distal_site/distal_site.csv
```

The caller therefore should not need to manually specify either the RFdiffusion input directory or the distal-site output path. Both can be derived from the pipeline run directory.

The main configurable argument for this stage is the peptide chain:

```text
--peptide-chain A
```

The result is a single CSV describing the distal-site selection for the RFdiffusion designs in that pipeline run.


## ProteinMPNN/RosettaRelax 

I am working on mpnn/relax scripts. First, I added --force flags so I can run the command on the same directory. Then, I added portability (environment/relative scripts). Exposed rounds as an arg. 

```sh
# the script is just a wrapper and will
# activate proteinmpnn/pyrosetta
conda activate biopython

python /home/dkouv/work/tfr1-binder-design/scripts/run_mpnn_relax_all.py \
  --pdb-dir /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/rfdiffusion \
  --scores /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/distal_site/distal_site.csv \
  --output-dir /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/mpnn_relax \
  --xml /home/dkouv/work/tfr1-binder-design/RosettaRelax/fast_relax_binder_A.xml \
  --proteinmpnn-dir /home/dkouv/ProteinMPNN \
  --force \
  --rounds 1
```

## AfCyc Oracle

Similar to ProteinMPNN/Relax, it was made more portable. There was an issue with my env w/r/t to GPU availability in JAX, also fixed. You can select which round to process, but it defaults to the latest round (usual case). We also need `XLA_PYTHON_CLIENT_PREALLOCATE=false` because JAX otherwise reserves too much of the GPU and immediately OOMs. 

```sh
cd /tmp
conda activate biopython

XLA_PYTHON_CLIENT_PREALLOCATE=false \
python /home/dkouv/work/tfr1-binder-design/scripts/run_afcyc_all.py \
  --design-dir /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/mpnn_relax \
  --output-dir /home/dkouv/work/tfr1-binder-design/experiments/003_workflow_dev_rfdiffusion/outputs/piecemeal-project/oracle \
  --target-pdb /home/dkouv/work/tfr1-binder-design/data/PDBs/6wrw_ChA_189-383.pdb \
  --params /home/dkouv/alphafold \
  --afcyc-env afcyc \
  --round 1 \
  --recycles 1
```

## Simple pipeline 

To get my head around all the ins and outs, we are going to run a simple shell script just to see if things chain together. It is 

