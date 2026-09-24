# 002: Apical Domain Beta Pairing

Set up the pyrosetta + torch env on WSL. 

The script to calculate the SS/ADJ tensors is in: `scripts/interface_tensors/make_interface_tensor.py`.

The primary goal of the experiment is if beta pairing causes the binder designs to fall apart somehow. The secondary goal is how to best design the tensor themselves.

I am not sure what the optimal settings for these are yet or how they can best be evaluated. Ultimately, we want our oracle to concur. 

However, that is very far downstream from backbone generation. Systematic evaluation of whether it was this step, the backbone generation, or the sequence design is really tough and would require a lot of compute. I think this generally one of reasons people pursue a large number of backbones and don't often report a lot tuning (AFAIK), but that's only a guess. I don't think RFpeptides commented on this, but reading more RFdiffusion papers would give me better context.  

My first assessment is this: 

1. If we want a diversity of backbones, we should try to limit the target adjacency to be the most prominent residues for backbone complementarity. The logic being that we let the model explore as much as it can while keeping only what we need. 

This is 210-212 in the target visually, but could be argued to 208-212. I will make a figure about this/copy the PSE here. AFAIK, the script has an off by one bug, I wrote about it before in the other repo. But this is why we have `--target_adj A211-213`.

2. The binder length is maybe better thought of as being more driven by wet lab experimental consideration. At least, it sets an upper bound. RFpeptides did 12-18-mers. Ideally, I would actually like to trial that range, the compute gets hard again however. It's very challenging to reason about if we'd like something smaller or larger - there are many nondecomposable effects on inter- and intra- molecular interactions. I think we may want to start towards 18 and see if the model is making the most of its residues, so to say. After inspecting those, I might have a better perspective on where to go next. We will take `--binderlen 18` as a first pass.

3. The amount of beta strand character required is probably 3, to pair with with the target's beta strand residues. Giving us `--binder_ss_len 4`.

This example script is relative to root, but we are in `experiments/002_beta_pair_apical`.

```sh
conda activate pyrosetta
# suspected an off by 1 error due to matrix assignment,
# so, --target_adj A211-213 instead of A210-212
python scripts/interface_tensors/make_interface_tensor.py \
  --input_pdb data/PDBs/6wrw_ChA_189-383.pdb \
  --out_dir experiments/002_beta_pair_apical/outputs/interface_tensors_A211-213_bl-18_bsl-3 \
  --binderlen 18 \
  --target_adj A211-213 \
  --binder_ss E \
  --binder_ss_len 3
```
Ran it from root. 

We have: `experiments/002_beta_pair_apical/outputs/interface_tensors_A211-213_bl-18_bsl-3`

Output appears okay. 

Now, we will start RFdiffusion. Too many args to make a sane dir name, so this is exp 1: `rfd_crop_tensor_1`

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv
cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_1" \
  inference.num_designs=10 \
  inference.input_pdb="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  'contigmap.contigs=[18-18 A189-383/0]' \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A211-213_bl-18_bsl-3" \
  scaffoldguided.scaffold_list="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A211-213_bl-18_bsl-3.txt" \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A210,A211,A212]'
```

~ 30 min.

Mistake in the script, moved to: `/home/dkouv/work/tfr1-binder-design/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_1`.

```sh
# from /Users/lkv206/work/tfr1-binder-design/experiments/002_beta_pair_apical
pymol outputs/rfd_crop_tensor_1/rfd_crop_tensor_1_*.pdb
```
Chain A: Binder
Chain B: Apical Domain TfR1 

Within PyMol, you can make them into states: 

```pml
python
from glob import glob
import re

files = sorted(
    glob("outputs/rfd_crop_tensor_1/rfd_crop_tensor_1_*.pdb"),
    key=lambda f: int(re.search(r'_(\d+)\.pdb$', f).group(1))
)

for f in files:
    cmd.load(f, "rfd_designs")
python end
```

View sticks: 

```sh
show sticks, rfd_designs and chain A
select target_near_binder, rfd_designs and chain B within 5 of (rfd_designs and chain A)
show sticks, target_near_binder
```

Scaled up to 100. 

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

mkdir -p "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_2"

conda activate SE3nv
cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_2/rfd_crop_tensor_2" \
  inference.num_designs=100 \
  inference.input_pdb="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  'contigmap.contigs=[18-18 A189-383/0]' \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A211-213_bl-18_bsl-3" \
  scaffoldguided.scaffold_list="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A211-213_bl-18_bsl-3.txt" \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A210,A211,A212]'
```

~ 3.5 hours. 

Many of these design look visually very odd compared to to what I had previously observed.  

I ran a script to investigate beta strand pairing. The env set up is here:  `experiments/002_beta_pair_apical/biotite_env_setup.md`.

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs outputs/rfd_crop_tensor_2 \
    --output-dir outputs/rfd_crop_tensor_2/beta_analysis \
    --target-chain B \
    --target-interface 208-212 \
    --binder-chain A
```

0% inter B pairing, 4% intra (100 designs). Script seems correct, but to check, I put previous beta pairs that looked good from previous campaign (https://github.com/DKchemistry/brain-peptide-shuttle-design) in `experiments/002_beta_pair_apical/inputs/15-mer-macrocycles`. 


```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs inputs/15-mer-macrocycles \
    --output-dir inputs/15-mer-macrocycles/beta_analysis \
    --target-chain B \
    --target-interface 208-212 \
    --binder-chain A
```
This output didn't seem right, 0% inter and 40% intra. Intra makes sense. Inter visibly looks better. It seems to be an issue from the *previous step*. The target PDB is renumbered. I will adjust the interface args, but to prevent more downstream problems, this should be fixed in a pipeline. 

In the cropped designs, it's at about 39-42.  

In the uncropped designs, it's at about 103-105, maybe 102-107 

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs inputs/15-mer-macrocycles \
    --output-dir inputs/15-mer-macrocycles/beta_analysis \
    --target-chain B \
    --target-interface 102-107 \
    --binder-chain A
```

```
Successfully analyzed: 10
Inter beta-sheet complementarity: 100.0% (10/10)
Intra beta-sheet complementarity: 40.0% (4/10)
Figure: inputs/15-mer-macrocycles/beta_analysis/beta_sheet_complementarity.png
```

That seems a little high. Certainly `tfr1_beta4_cyc_chain_a_2.pdb` looks terrible in PyMol.

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs inputs/15-mer-macrocycles \
    --output-dir inputs/15-mer-macrocycles/beta_analysis \
    --target-chain B \
    --target-interface 103-105 \
    --binder-chain A
```

Successfully analyzed: 10
Inter beta-sheet complementarity: 100.0% (10/10)
Intra beta-sheet complementarity: 40.0% (4/10)
Figure: inputs/15-mer-macrocycles/beta_analysis/beta_sheet_complementarity.png

I feel it should be 9/10 at best. PyMol does not even draw a beta strand in `tfr1_beta4_cyc_chain_a_2.pdb`'s binder. It's drawn as a floppy loop. 

The default minimal pairing is not exposed as an arg, but it is 2 in the above. I changed it to 3 for this run below. I will update the pairing as an arg to avoid confusion later: 

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs inputs/15-mer-macrocycles \
    --output-dir inputs/15-mer-macrocycles/beta_analysis \
    --target-chain B \
    --target-interface 103-105 \
    --binder-chain A
```

Successfully analyzed: 10
Inter beta-sheet complementarity: 80.0% (8/10)
Intra beta-sheet complementarity: 40.0% (4/10)
Figure: inputs/15-mer-macrocycles/beta_analysis/beta_sheet_complementarity.png

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs inputs/15-mer-macrocycles \
    --output-dir inputs/15-mer-macrocycles/beta_analysis \
    --target-chain B \
    --target-interface 103-105 \
    --binder-chain A \
    --min-inter-beta-pairs 3 \
    --min-intra-beta-pairs 2 
```

Successfully analyzed: 10
Inter beta-sheet complementarity: 80.0% (8/10)
Intra beta-sheet complementarity: 40.0% (4/10)
Figure: inputs/15-mer-macrocycles/beta_analysis/beta_sheet_complementarity.png

Note: Forgot, they are 14 member not 15 member. 

Let's return to the current designs. We will start very permissive on the target interface and min inter pairs to see if the b pairing worked at all.

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs outputs/rfd_crop_tensor_2 \
    --output-dir outputs/rfd_crop_tensor_2/beta_analysis \
    --target-chain B \
    --target-interface 38-43 \
    --binder-chain A \
    --min-inter-beta-pairs 2 \
    --min-intra-beta-pairs 2 
```

Successfully analyzed: 100
Inter beta-sheet complementarity: 91.0% (91/100)
Intra beta-sheet complementarity: 4.0% (4/100)
Figure: outputs/rfd_crop_tensor_2/beta_analysis/beta_sheet_complementarity.png

More strict on inter requirements:

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs outputs/rfd_crop_tensor_2 \
    --output-dir outputs/rfd_crop_tensor_2/beta_analysis \
    --target-chain B \
    --target-interface 38-43 \
    --binder-chain A \
    --min-inter-beta-pairs 3 \
    --min-intra-beta-pairs 2 
```

Successfully analyzed: 100
Inter beta-sheet complementarity: 50.0% (50/100)
Intra beta-sheet complementarity: 4.0% (4/100)
Figure: outputs/rfd_crop_tensor_2/beta_analysis/beta_sheet_complementarity.png


I would argue this is even too kind based on visual inspection, but some of that is clouded by looking at PyMol to judge. Regardless, it's clearly worse than the previous 14-mer run. 

I think we need to set up new tensors, RFd 14-mers again at larger scale, set up the interface residues, and see what is going on. You can argue 6WRW has four b-pairs, but they are not contigious. I'll keep `--binder_ss_len 3` but this could be changed. This idea is that RFdiffusion will find a good binding mode that may be >3, but should be at least 3. `--target_adj` should probably be expanded, as low as 208 is B-paired in 6WRW imo. So 209 due to the script issue. 

Run from root:

```sh
conda activate pyrosetta
# suspected an off by 1 error due to matrix assignment,
# so, --target_adj A209-213 instead of A208-212
python scripts/interface_tensors/make_interface_tensor.py \
  --input_pdb data/PDBs/6wrw_ChA_189-383.pdb \
  --out_dir experiments/002_beta_pair_apical/outputs/interface_tensors_A209-213_bl-14_bsl-3 \
  --binderlen 14 \
  --target_adj A209-213 \
  --binder_ss E \
  --binder_ss_len 3
```

# `rfd_crop_tensor_3` - 10 designs

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv
cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_3/rfd_crop_tensor_3" \
  inference.num_designs=10 \
  inference.input_pdb="$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  'contigmap.contigs=[14-14 A189-383/0]' \
  inference.model_runner=ScaffoldedSampler \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.scaffold_dir="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A209-213_bl-14_bsl-3" \
  scaffoldguided.scaffold_list="$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A209-213_bl-14_bsl-3.txt" \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A208,A209,A210,A211,A212]'
```

```sh
conda activate biotite
python ../../scripts/beta_sheet_complementarity.py \
    --input-pdbs outputs/rfd_crop_tensor_3 \
    --output-dir outputs/rfd_crop_tensor_3/beta_analysis \
    --target-chain B \
    --target-interface 38-43 \
    --binder-chain A \
    --min-inter-beta-pairs 2 \
    --min-intra-beta-pairs 2 
```

Successfully analyzed: 10
Inter beta-sheet complementarity: 0.0% (0/10)
Intra beta-sheet complementarity: 30.0% (3/10)
Figure: outputs/rfd_crop_tensor_3/beta_analysis/beta_sheet_complementarity.png

I am not sure what to make of that. Visually, PyMol isn't assigning the beta sheet cartoon, but some of these backbones do appear to have inter sheet complementarity. Have I made some mistake in my workflow? Should `'contigmap.contigs=[14-14 A189-383/0]'` and `'ppi.hotspot_res=[A208,A209,A210,A211,A212]'` have been changed b/c of PDB renumbering? Maybe, but the binders **are** at the right spot. RFdiffusion has a propensity for alpha helicies and we are not seeing them, so the SS/ADJ method is too some extent working. We are even getting some visually "non-awful" B-pairs. I don't understand how the designs from the last project were so much better (what I've called the 15-mer set). That workflow is here: https://github.com/DKchemistry/brain-peptide-shuttle-design/blob/main/METHODS.md

I will match it as exactly as I can, I just can't do full length chain A due to VRAM. 

# `rfd_crop_tensor_4` - `rfd_crop_tensor_3` had a mistake! 

There was a mistake in the downstream analysis of `rfd_crop_tensor_3`: `--target-interface 38-43` used the wrong RFdiffusion output numbering. For this experiment, the conditioned TfR1 residues are original A209–A212, which become **B35–B38** in the RFdiffusion output.

## Generate the tensors

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate pyrosetta

python "$ROOT/scripts/interface_tensors/make_interface_tensor.py" \
  --input_pdb "$ROOT/data/PDBs/6wrw_ChA_189-383.pdb" \
  --out_dir "$ROOT/experiments/002_beta_pair_apical/outputs/interface_tensors_A210-213_bl-14_bsl-4" \
  --binderlen 14 \
  --target_adj A210-213 \
  --binder_ss E \
  --binder_ss_len 4
```

`A210-213` is intentional: because of the off-by-one bug in `make_interface_tensor.py`, this conditions actual target residues **A209–A212**.

## Run RFdiffusion

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
RFD=/home/dkouv/RFdiffusion

conda activate SE3nv

mkdir -p "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4"

cd "$RFD"

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix="$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4/rfd_crop_tensor_4" \
  inference.num_designs=10 \
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

## Analyze beta-sheet complementarity

For these outputs, original TfR1 residues **A209–A212** are renumbered to **B35–B38**.

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate biotite

python "$ROOT/scripts/beta_sheet_complementarity.py" \
  --input-pdbs "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4" \
  --output-dir "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4/beta_analysis" \
  --target-chain B \
  --target-interface 35-38 \
  --binder-chain A \
  --min-inter-beta-pairs 2 \
  --min-intra-beta-pairs 2
```
Successfully analyzed: 10
Inter beta-sheet complementarity: 90.0% (9/10)
Intra beta-sheet complementarity: 30.0% (3/10)
Figure: /home/dkouv/work/tfr1-binder-design/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4/beta_analysis/beta_sheet_complementarity.png

Okay, much more sane. 

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate biotite

python "$ROOT/scripts/beta_sheet_complementarity.py" \
  --input-pdbs "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4" \
  --output-dir "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4/beta_analysis" \
  --target-chain B \
  --target-interface 35-38 \
  --binder-chain A \
  --min-inter-beta-pairs 3 \
  --min-intra-beta-pairs 2
```

Successfully analyzed: 10
Inter beta-sheet complementarity: 50.0% (5/10)
Intra beta-sheet complementarity: 30.0% (3/10)
Figure: /home/dkouv/work/tfr1-binder-design/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_4/beta_analysis/beta_sheet_complementarity.png

Also reasonable. 

Visual inspection of the designs looks reasonable. What is hard to tell presently is the relationship of `'ppi.hotspot_res=[A209,A210,A211,A212]'` to the parameters by which we set the SS/ADJ tensors? They both, in a sense, condition a hot spot or a hot spot-like behavior. For example, there is a neighboring alpha helix to the apical beta strand that could plausibly be used to stabilize interactions of the macrocycle, spanning roughly 166 - 173 in the renumbered PDB. Is it possible to hot spot at the helix while Beta pairing at the strand? Does the internals of RFdiffusion even theoretically support that? 

It seems too. 

But, to clarify, the mistake here has been the analysis (at least in part) due to renumbering: 

1. Cropped input PDB
   A189–A383

2. Tensor/RFdiffusion input arguments
   still use A189–A383 numbering
   (plus the tensor script's separate +1 bug)

3. RFdiffusion output
   binder A1–A14
   target B15–...

So let's re run: `rfd_crop_tensor_2`.

This should have been about 38-43.

```sh
ROOT=/home/dkouv/work/tfr1-binder-design
conda activate biotite
python "$ROOT/scripts/beta_sheet_complementarity.py"  \
    --input-pdbs $ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_2 \
    --output-dir $ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_2/beta_analysis \
    --target-chain B \
    --target-interface 38-43 \
    --binder-chain A \
    --min-inter-beta-pairs 3 \
    --min-intra-beta-pairs 2
```

Successfully analyzed: 100
Inter beta-sheet complementarity: 50.0% (50/100)
Intra beta-sheet complementarity: 4.0% (4/100)
Figure: /home/dkouv/work/tfr1-binder-design/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_2/beta_analysis/beta_sheet_complementarity.png


So it's an estimate presently, but I don't think macrocycle length has a huge impact on inter-strand. I am guessing it does on intra. 

So I think there are a few design strategies: 

14-mer -> try to stack two intra-strands with the inter-strand. 
* It may be difficult to get to the helix at this macrocycle length/geometry. 
* This is essentially `rfd_crop_tensor_4`. I want to scale this to a 100. 
* This is`rfd_crop_tensor_5` 

14-mer -> try to only stack the inter-strand and reach for the helix via hot spots. 
* It maybe tough but it could work? Won't know until we try. 
* could be `rfd_crop_tensor_6`

18-mer -> try to only stack the inter-strand and reach for the helix via hot spots.
* Seems to be the opposite case of the 14-mer, intra will be (my hypothesis) tough because diffusion has more options to spread out. 
* inter may be easier by the same logic. 
* could be `rfd_crop_tensor_7`

18-mer -> try to stack the intra strands with the inter-stand, 
* This is basically to stress test my hypothesis. 
* could be `rfd_crop_tensor_7`

While these are all good ideas, they are putting the cart before the horse. I think it is better to get a pipeline set up through the rest of the design process. Later, it would be easier to intervene in the pipeline to try other strategies and compare what I'd like to at the oracle itself.

# `rfd_crop_tensor_5`

Visually, the designs look reasonable in PyMol.

```sh
ROOT=/home/dkouv/work/tfr1-binder-design

conda activate biotite

python "$ROOT/scripts/beta_sheet_complementarity.py" \
  --input-pdbs "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_5" \
  --output-dir "$ROOT/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_5/beta_analysis" \
  --target-chain B \
  --target-interface 35-38 \
  --binder-chain A \
  --min-inter-beta-pairs 3 \
  --min-intra-beta-pairs 2
```
Successfully analyzed: 100
Inter beta-sheet complementarity: 54.0% (54/100)
Intra beta-sheet complementarity: 10.0% (10/100)
Figure: /home/dkouv/work/tfr1-binder-design/experiments/002_beta_pair_apical/outputs/rfd_crop_tensor_5/beta_analysis/beta_sheet_complementarity.png

Seems reasonable as well. I am a little surprised the intra-strand complementarity is this low. 