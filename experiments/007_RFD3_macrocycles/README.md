# Background:

This experiment concerns my initial implementation of cyclic positional encoding in RFD3 and testing generation capabilities. 

Repo Link: https://github.com/DKchemistry/foundry/tree/feat/rfd3-cyclic-positional-encoding

Commit SHA: 4bdb42af496f7e77c92a2db85afeb77788ea204b

## Install: 

On my WSL2 system:

```sh
git clone https://github.com/DKchemistry/foundry.git
cd foundry

git fetch origin
git switch feat/rfd3-cyclic-positional-encoding
```

Pin to the inital implementation commit. 

```sh
git checkout 4bdb42af496f7e77c92a2db85afeb77788ea204b
```

Use uv: 

```sh
uv python install 3.12
uv venv --python 3.12
source .venv/bin/activate
uv pip install -e '.[all,dev]'
```
This implementation will be tested with Python 3.12, Torch 2.14.0+cu130, CUDA working on GTX 1660 Ti.

Install checkpoint: 

```sh
foundry install rfd3
```

This was already satisfied.

## Update: Oct 4th, 2026 ~7:30 PM CPH Time 

I updated my implementation of RFD3 to use unique residues to calculate cyclic chain length rather than just counting along the index to match RF3 more closely. This should not affect the downstream results but re-running the analysis is not a bad idea. 

Initial implementation:
4bdb42af496f7e77c92a2db85afeb77788ea204b

Updated CRPE implementation:
69c70a88e0ed7348602ecfcef55d9ee6082cd675

## Test 1: Default monomer and binder design

### Binder/PPI

First, I will assess if my changes affected default behavior in monomer and binder design. 

I will use the [PPI example](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/examples/protein_binder_design.md) given in the tutorial. 


```sh
# pwd = $HOME/foundry
python "$HOME/foundry/models/rfd3/src/rfd3/run_inference.py" \
  out_dir="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/test_1/binder" \
  ckpt_path="$HOME/.foundry/checkpoints/rfd3_latest.ckpt" \
  inputs="$HOME/foundry/models/rfd3/docs/examples/protein_binder_design.json" \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2
```
Log appears normal relative to other runs I have done with RFD3. One potential concern: 

`WARNING:rfd3.model.layers.layer_utils:[rank: 0] Using nn.RMSNorm instead of apex.normalization.fused_layer_norm.FusedRMSNorm.Ensure you're using the correct apptainer`

But I am not sure I introduced this. I will check later. 

Run time is normal on my WSL install, I have run this test previously when installing RFD3.

Run time: 32m 15s

Output is normal by visual inspection in PyMol. Structures appear like reasonable binders. No accidental cyclicity was introduced. 

### Monomer

I made my own monomer example. See: `input_jsons/monomer_design.json`. It requests 20-40 AA monomers and otherwise uses sane defaults. 

```sh
# pwd = $HOME/foundry
python "$HOME/foundry/models/rfd3/src/rfd3/run_inference.py" \
  out_dir="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/test_1/monomer" \
  ckpt_path="$HOME/.foundry/checkpoints/rfd3_latest.ckpt" \
  inputs="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/input_jsons/monomer_design.json" \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2
```

Run time: 1m 20s

Output is normal by visual inspection in PyMol. Structures appear like reasonable monomers. No accidental cyclicity was introduced. 

## Test 2: Default macrocycle monomer and binder design.

### Monomer 

I am following what was reported in RFpeptides, but scaling down due to compute budget: 

> We added the cyclic positional encoding scheme to RFdiffusion and observed robust generation of diverse macrocyclic peptides (Fig. 1b,c and Supplementary Fig. 2). Similar to the previously described work on designing monomeric cyclic peptides with physics-based methods7 and AfCycDesign11, we observed 9,045 and 8,913 structurally unique 10-residue and 12-residue backbones, respectively, when 48,000 macrocycle backbones were generated for each size (Supplementary Fig. 2). The distribution of phi and psi values in these generated backbones is similar to the standard Ramachandran plot for protein structures (Supplementary Fig. 2), suggesting that generated backbones do not require extensive d-amino acids to stabilize the generated structures7. While we did not attempt to comprehensively enumerate the structural space of cyclic peptide monomers, RFpeptides can readily be scaled up to comprehensively cover the structural space accessible to macrocyclic peptides. 

It can be seen from the figures that folding the inverse design sequences also occured. We already have AFCyc as an oracle on this machine. I first need to see what our throughput and storage requirements are. I will start by making 8 monomers of 10-mers and 12-mers. 

/home/dkouv/foundry/models/rfd3/docs/examples/macrocycle_monomer.json


```sh
# pwd = $HOME/foundry
python "$HOME/foundry/models/rfd3/src/rfd3/run_inference.py" \
  out_dir="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/test_2/macrocycle_monomer" \
  ckpt_path="$HOME/.foundry/checkpoints/rfd3_latest.ckpt" \
  inputs="$HOME/foundry/models/rfd3/docs/examples/macrocycle_monomer.json" \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2
```

1m 33s

The macrocycles look visually correct in PyMol, there is N -> C connectivity. 

Let's make 500 of each by setting `diffusion_batch_size` to 500. I am not sure what the difference between that and increasing `n_batches` is. Both lead to more designs overall, but whether they mean scientifically different things I am not sure. I am also not sure about `step_scale` and `gamma`, but these were recommended in binder design. I will ask in the discord or read more.

```sh
python "$HOME/foundry/models/rfd3/src/rfd3/run_inference.py" \
  out_dir="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/test_2/macrocycle_monomer_500" \
  ckpt_path="$HOME/.foundry/checkpoints/rfd3_latest.ckpt" \
  inputs="$HOME/foundry/models/rfd3/docs/examples/macrocycle_monomer.json" \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2 \
  diffusion_batch_size=500
```

15m 3s

RFD3 is very fast. I will make 10,000 of each and pause there. That should be enough to get some meaningful signal. It may even be possible to get to 50K (est. 25 hr of compute), but I am not comfortable running that on the laptop. 10,000 should be 5 hours of compute, which is not small. The sampler settings are debateable. The `step_scale=3` and `gamma_0=0.2` are recommended for the designability of binders, but comes at some loss of diversity. The base settings of RFD3 (`step_scale=1.5` and `gamma_0=0.6`) encourage more diversity. I think establishing the baseline is a good initial experiment. I will also use `is_non_loopy=True`. I was likely far too optimistic about RFD3 throughput, though I hope I was wrong. 


```sh
python "$HOME/foundry/models/rfd3/src/rfd3/run_inference.py" \
  out_dir="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/test_2/macrocycle_monomer_10k" \
  ckpt_path="$HOME/.foundry/checkpoints/rfd3_latest.ckpt" \
  inputs="$HOME/work/tfr1-binder-design/experiments/007_RFD3_macrocycles/test_2/macrocycle_monomer_10K/macrocycle_monomer_10K.json" \
  diffusion_batch_size=500 \
  n_batches=20
```
4h 52m 18s

This is about 50% VRAM but already 98% GPU util, so further increments of `diffusion_batch_size` may not be offering much more gain.

I can scale more later if required, but I think it is better to try to reproduce the analysis done in the RFpeptides work.

### Monomer: Cluster Diversity and Phi/Psi Angles as a Function of Designs 

In RFpeptides, one of the methods used to assess the diversity of design space is by Calpha RMSD clustering. If cluster size grows as a function of designs, we can expect that the method is generating a diversity of Calpha designs. I perform a similar experiment in `experiments/007_RFD3_macrocycles/test_2/macrocycle_monomer_10k/clustering`, but instead using C-backbone RMSD. 

The basic approach is to use cyclic permutations to assign a psuedo-"low energy bin" cluster. Consider a five membered cyclic peptide composed of A1 through A5 and a second five membered cyclic peptide composed of B1 through B5. We perform the superposition and calculate RMSD of A1-A5 onto B1-B2-B3-B4-B5, B2-B3-B4-B5-B1, B3-B4-B5-B1-B2, etc etc. Imagine one permutation of B produces an RMSD < 0.5A, A and B then become members of the same cluster. We try A1-A5 against permutations of C1-C5, we do not need to try B onto C as it already belongs in the first cluster with A. C may not have a permutation that meets our threshold, so instead we begin using C to define cluster 2 and align it to D, this process repeats until we have evaluated every potential cluster representative. If designs are sampling diverse C-backbone RMSDs, we expect clusters to increase monotonically with design size. That is indeed what we see: 

![alt text](test_2/macrocycle_monomer_10k/clustering/figures/macrocycle_cluster_growth.png)

Relevant code: `experiments/007_RFD3_macrocycles/scripts/cluster_macrocycle_backbones.py`
Relevant code: `experiments/007_RFD3_macrocycles/scripts/plot_macrocycle_clustering.py`

That does not, in and of itself, guarantee anything about the reasonability of the backbones. So, our next experiment is to calculate the Phi/Psi angle distributions of our designs to assess whether we are sampling reasonable conformational space. This is, to me, still mostly a qualitative check. 

![alt text](test_2/macrocycle_monomer_10k/ramachandran/figures/macrocycle_ramachandran.png)

Relevant Code: `experiments/007_RFD3_macrocycles/scripts/calculate_macrocycle_torsions.py`
Relevant Code: `experiments/007_RFD3_macrocycles/scripts/plot_macrocycle_ramachandran.py`