# Cropping TfR1

I watched some RFdiffusion tutorials from the Rosetta team. It is advised to crop your target protein in binder design to improve throughput. This is treated as more of an art than a science it seems. 

I am using TfR1 from PDB 6WRW. 

It is in: `data/PDBs`. I am in: `experiments/001_crop_tfr1`.

## Crop to the apical domain

I cropped residues 189-383 of Chain A in 6wrw via PyMol (194 res total). Saved here: `experiments/001_crop_tfr1/input/6wrw_ChA_189-383.pdb`. Reopening it in PyMol suggests all the PDB metadata is as expected: starts at res 189, ends at res 383, still Chain A. 

## RFdiffusion

I will mimic the macrocycle binder design I will want later. 

This code will run on my WSL machine with a 1660 Ti. I've profiled RFdiffusion before at roughly this residue number (~200 ish): 

| Target | Time | RFdiffusion VRAM | Total peak VRAM |
|---:|---:|---:|---:|
| 50 residues | 1 min 26 sec | 1.92 GiB | 2.39 GiB |
| 100 residues | 1 min 27 sec | 2.29 GiB | 2.73 GiB |
| 200 residues | 3 min 14 sec | 2.67 GiB | 3.10 GiB |
| 300 residues | 6 min 11 sec | 4.70 GiB | 5.13 GiB |
| 400 residues | 10 min 43 sec | 4.73 GiB | 5.16 GiB |
| 450 residues | 13 min 8 sec | 5.13 GiB | 5.51 GiB |
| 500 residues | 15 min 38 sec | 5.36 GiB | 5.73 GiB |

So it should hopefully execute fine. 

I will try this first. I don't recall if the CUDA flag needs to be explicit. This should be the correct `contigmap.contigs` and `ppi.hotspot_res` syntax I have used before. The directory should be synced using rsync. 

```sh
conda activate SE3nv
cd /home/dkouv/RFdiffusion

time python scripts/run_inference.py \
  --config-name base \
  inference.output_prefix=/home/dkouv/work/tfr1-binder-design/experiments/001_crop_tfr1/output/6wrw_ChA_189-383 \
  inference.num_designs=1 \
  inference.input_pdb=/home/dkouv/work/tfr1-binder-design/experiments/001_crop_tfr1/input/6wrw_ChA_189-383.pdb \
  'contigmap.contigs=[15-15 A189-383/0]' \
  inference.cyclic=True \
  inference.cyc_chains='a' \
  diffuser.T=50 \
  'ppi.hotspot_res=[A209,A210,A211,A212]'
```
~ 3 minutes to execute.

I aligned them in PyMol. The RMSD is tiny, 0.09. The binder diffused to the right hot spots. 

```sh
align 6wrw_ChA_189-383, 6wrw_ChA_189-383_0
 Match: read scoring matrix.
 Match: assigning 195 x 210 pairwise scores.
 MatchAlign: aligning residues (195 vs 210)...
 MatchAlign: score 1009.000
 ExecutiveAlign: 780 atoms aligned.
 ExecutiveRMS: 46 atoms rejected during cycle 1 (RMSD=0.19).
 ExecutiveRMS: 41 atoms rejected during cycle 2 (RMSD=0.13).
 ExecutiveRMS: 34 atoms rejected during cycle 3 (RMSD=0.11).
 ExecutiveRMS: 17 atoms rejected during cycle 4 (RMSD=0.10).
 ExecutiveRMS: 8 atoms rejected during cycle 5 (RMSD=0.09).
 Executive: RMSD =    0.090 (634 to 634 atoms)
```

It seems we can safely crop to this ~ 200 residue system without much worry. 
