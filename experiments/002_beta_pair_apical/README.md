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
