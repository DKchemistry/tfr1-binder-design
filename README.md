# TFR1 Binder Design

## Experiments 

### 001: Cropping TfR1 to the Apical Domain

`experiments/001_crop_tfr1`

RFdiffusion scales agressively with residue count. This experiment looks at how we can crop TfR1 to improve throughput without corrupting the structure.

In a single run of RFdiffusion with only the apical domain residues of TfR1 as input (res. 189-383) as a target and 15-mer macrocycle as binder, we recover the target structure with low RMSD (0.09) and with the intended hot spot interactions. Roughly ~3 min to execute. Night and day improvement over 641 residue/Apple MPS runs (~2 hrs). 

### 002: Apical Domain Beta Pairing

`experiments/002_beta_pair_apical`

I previously used a β-pair targeting conditioning method to encourage sampling of β-strand interactions between TfR1/macrocycle ([Sappington et al., 2026](https://doi.org/10.1038/s41467-025-67866-3)). This seemed successful for full length Chain of 6WRW. The scripts needed to calculate the SS/ADJ tensor blocks were downloaded here: `scripts/interface_tensors/`. This experiment examines if this behavior holds at the apical domain crop of TfR1.



