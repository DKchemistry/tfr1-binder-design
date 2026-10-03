# Macrocycle backbone torsion analysis

This directory contains the phi, psi, omega, peptide-bond, and Ramachandran
validation analysis for the 10,000 generated 10-mers and 10,000 generated
12-mers. It is intentionally separate from the experiment README.

## Scope

The RFpeptides-like primary figure pools every phi/psi observation for a given
peptide length. It does not separate structures by RMSD cluster, matching the
presentation in RFpeptides Supplementary Figure S2I.

Every residue in each cyclic peptide contributes one complete set of torsions:

- phi: `C(previous)-N-CA-C`
- psi: `N-CA-C-N(next)`
- omega: `CA-C-N(next)-CA(next)`

The residue before residue 1 is the final residue, and the residue after the
final residue is residue 1. The closure bond is therefore evaluated explicitly
rather than being left undefined by a linear-chain helper.

## Software and reference data

The structures are parsed with Biotite 1.6.0 and the torsions are calculated
with `biotite.structure.dihedral`. Biotite is a reasonable substitute for
Biopython here and is preferred because it is already familiar to the project
owner. The only important implementation detail is that
`dihedral_backbone()` is a linear-chain convenience function, so this analysis
uses the general four-coordinate `dihedral()` function to include the cyclic
wraparound.

Ramachandran validation uses the six MolProbity Top8000 residue classes:

1. general
2. glycine
3. cis-proline
4. trans-proline
5. pre-proline
6. isoleucine or valine

Proline classification uses the omega angle of the peptide bond entering the
proline. The locally stored Top8000 grids are unchanged copies of the official
Richardson Lab `rotarama_data` repository at commit
`76aae74c6e1f834775f8df4700c79602c6a8b9ec`. See
`reference_data/molprobity_top8000/PROVENANCE.md` for checksums and licensing.

Top8000 scores use bilinear interpolation between 2-degree bin centers with
periodic wrapping. Classification uses the current cctbx `ramalyze` cutoffs:

- favored: score greater than or equal to 0.02 for every class
- allowed general: score greater than or equal to 0.0005
- allowed cis-proline: score greater than or equal to 0.002
- allowed all other classes: score greater than or equal to 0.001
- outlier: below the class-specific allowed cutoff

## Commands

Run in the Biotite conda environment from the repository root:

```sh
MPLCONFIGDIR=/tmp/rama-matplotlib \
  /home/dkouv/miniforge3/envs/biotite/bin/python \
  experiments/007_RFD3_macrocycles/scripts/calculate_macrocycle_torsions.py \
  --overwrite

MPLCONFIGDIR=/tmp/rama-matplotlib \
  /home/dkouv/miniforge3/envs/biotite/bin/python \
  experiments/007_RFD3_macrocycles/scripts/validate_macrocycle_ramachandran.py \
  --overwrite

MPLCONFIGDIR=/tmp/rama-matplotlib \
  /home/dkouv/miniforge3/envs/biotite/bin/python \
  experiments/007_RFD3_macrocycles/scripts/plot_macrocycle_ramachandran.py

MPLCONFIGDIR=/tmp/rama-matplotlib \
  /home/dkouv/miniforge3/envs/biotite/bin/python \
  experiments/007_RFD3_macrocycles/scripts/plot_macrocycle_omega.py
```

## Validation performed on 2026-10-03

### Dataset completeness and parsing

- 20,000 structures were discovered: 10,000 10-mers and 10,000 12-mers.
- Exactly 100,000 10-mer and 120,000 12-mer residue records were written.
- Every structure contained one amino-acid chain with the residue count implied
  by its filename.
- Every residue contained exactly one `N`, `CA`, and `C` atom.
- All 660,000 calculated phi, psi, and omega angles were finite.

### Independent torsion cross-check

A deterministic sample of 100 structures (50 of each length) was reparsed with
Biopython 1.88. All 3,300 phi/psi/omega values were recalculated independently
with `Bio.PDB.vectors.calc_dihedral`.

- maximum circular difference from Biotite: 0.000026530 degrees
- mean circular difference: 0.000007043 degrees
- every difference was below 0.0001 degrees

This confirms the atom ordering, angle convention, and cyclic wraparound
implementation independently of Biotite.

### Top8000 interpolation cross-check

The SciPy interpolation used by the validation script was compared with a
literal implementation of cctbx `NDimTable.valueAt()` at 10,000 random
phi/psi pairs.

- maximum absolute score difference: `2.22e-16`
- all values matched within an absolute tolerance of `1e-14`

### Ramachandran results

| Length | Favored | Allowed | Outlier | Structures with at least one outlier |
|---:|---:|---:|---:|---:|
| 10 | 99,491 (99.491%) | 463 (0.463%) | 46 (0.046%) | 42 (0.42%) |
| 12 | 119,692 (99.743%) | 293 (0.244%) | 15 (0.0125%) | 14 (0.14%) |

These values comfortably satisfy the conventional MolProbity targets of more
than 98% favored and fewer than 0.2% outliers. This is strong evidence that the
generated phi/psi values occupy established protein-like conformational
regions. It is not, by itself, a proof of sequence-specific foldability or
macrocycle stability.

### Peptide-bond and closure diagnostics

Mean C-to-next-N distances were 1.3203 A for 10-mers and 1.3172 A for 12-mers.
Mean closure-bond distances were 1.3238 A and 1.3220 A, respectively.

A deliberately broad 1.2-1.5 A diagnostic interval was also applied. This is
a conspicuous-geometry screen, not a proper residue-specific restraint
z-score.

| Length | All bonds outside screen | Closure bonds outside screen | Structures affected |
|---:|---:|---:|---:|
| 10 | 80 / 100,000 (0.080%) | 9 / 10,000 (0.090%) | 48 / 10,000 (0.48%) |
| 12 | 7 / 120,000 (0.0058%) | 1 / 10,000 (0.010%) | 6 / 10,000 (0.06%) |

Most of the conspicuous distances are therefore rare, and they are not
specific to the cyclic closure. A small number of 10-mer structures contain
several distorted bonds and should be treated as local geometry failures even
though the population-level Ramachandran result is excellent. The most extreme
example is `macrocycle_monomer_10K_macrocycle_10_5_model_114.cif.gz`, which has
six bonds outside the broad interval.

### Omega results

Using MolProbity's 30-degree cis/trans convention:

| Length | cis X-Pro | cis non-Pro | twisted X-Pro | twisted non-Pro |
|---:|---:|---:|---:|---:|
| 10 | 47 | 0 | 0 | 1 |
| 12 | 48 | 0 | 0 | 1 |

The two twisted non-proline bonds are individually identifiable in the torsion
table. The omega histogram in this directory uses all generated structures and
is labelled in peptide-bond observations. It is not labelled as clusters.

RFpeptides Supplementary Figure S2J instead appears to use cluster
representatives. Reproducing that panel faithfully requires a small extension
to the existing PyRosetta clustering workflow so that it permanently saves one
representative per cluster. That work should remain separate from the present
all-design torsion extraction so both scripts stay easy to understand.

## Output guide

- `all_structures_torsions.csv`: all cyclic torsions and C-N distances; large,
  reproducible, and intentionally ignored by git
- `validation_summary.json`: parsing, topology, distance, and omega summary
- `ramachandran_validation.csv`: per-residue Top8000 scores; large,
  reproducible, and intentionally ignored by git
- `ramachandran_quality.csv`: compact counts by length and residue class
- `ramachandran_validation_summary.json`: compact Top8000 provenance and results
- `figures/macrocycle_ramachandran.png`: RFpeptides S2I-style pooled plot
- `figures/macrocycle_ramachandran_by_residue_class.png`: diagnostic plot by
  MolProbity residue class
- `figures/macrocycle_omega.png`: all-design omega histogram on linear and log
  scales
