# TFR1 Binder Design

## Experiments 

### 001: Cropping TfR1 to the Apical Domain

`experiments/001_crop_tfr1`

RFdiffusion scales agressively with residue count. This experiment looks at how we can crop TfR1 to improve throughput without corrupting the structure.

Outcome: In a single run of RFdiffusion with only the apical domain residues of TfR1 as input (res. 189-383) as a target and 15-mer macrocycle as binder, we recover the target structure with low RMSD (0.09) and with the intended hot spot interactions. Roughly ~3 min to execute. Night and day improvement over 641 residue/Apple MPS runs (~2 hrs). 

### 002: Apical Domain Beta Pairing

`experiments/002_beta_pair_apical`

I previously used a β-pair targeting conditioning method to encourage sampling of β-strand interactions between TfR1/macrocycle ([Sappington et al., 2026](https://doi.org/10.1038/s41467-025-67866-3)). This seemed successful for full length Chain of 6WRW. The scripts needed to calculate the SS/ADJ tensor blocks were downloaded here: `scripts/interface_tensors/`. This experiment examines if this behavior holds at the apical domain crop of TfR1.

Outcome: Yes, the behavior remains. It seems to work as the paper implies at scales up to 100 backbones. This also highlighted an issue to watch out for during later workflows/analyses of RFdiffusion. When experimenting with 14-mer and 18-mer systems, I was getting odd results. It turns out that, in the diffusion output, when the binder is Chain A (which seems to be convention) and the target is Chain B, the Chain B will start at n+1 of the last residue in Chain A. So, in 14-mers, residue 15 of Chain B is the first residue. This can be very confusing in later analyses. For example, I was looking at β-pair complementarity after diffusion but my hotspot residues are not identical between 14-mer and 18-mer systems, which led me to misinterpert the results. I think this will become an issue later as well if I don't note it. 

### 003: Workflow development for RFdiffusion 

`experiments/003_workflow_dev_rfdiffusion`

During experiment 2, it became obvious that comparing any two conditions for some endpoint will require a computational pipeline with the ability to select relevant parameters at run time. This experiment focuses on developing such a pipeline. Many of the scripts from my previous project have been updated to be more portable and expose more arguments to the user. See `scripts/`. I have noted some issues on github regarding what is incomplete and potential footguns. The basic pipeline with all the scripts is being tested at the moment and will be revised over time as required. 

### 004: Oracle Success Rates: Scaling Beta Pairing 

To my knowledge, the combination of cyclic postional encoding (e.g. RFpeptides) and β-pair conditioning has not been previously reported. It would be interesting to set an initial baseline regarding whether this combination behaves as desired and, ideally, is supported by the oracle. This is challenging as my compute is still quite constrained. In the β-pairing work, *in silico* success rates varied from as low as 0.6% for FCRL5 and as high as 19.8% for α-CTX. I don't have a baseline for cropped TfR1 and would like to establish a baseline to later explore avenues to increase my success rates and my throughput. 
  





