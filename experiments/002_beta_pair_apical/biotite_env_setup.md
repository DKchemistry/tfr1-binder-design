Yes. I’d make this reproducible with **two small files**: a pinned `environment.yml` and a one-time setup script for the DSSP chemical-component dictionary.

The issue was twofold:

* `dssp 4.6.1` + current `libmcfp` gave an ABI/linker failure, so pinning DSSP to `4.5.6` avoided that.
* DSSP also needs the wwPDB `components.cif` dictionary, which conda did not populate automatically in this environment.

Use this `environment.yml`:

```yaml
name: biotite

channels:
  - conda-forge

dependencies:
  - python=3.11
  - biotite=1.6.0
  - dssp=4.5.6
  - matplotlib=3.10
  - numpy=2.3
  - curl
```

Then a setup script, e.g. `setup_biotite_env.sh`:

```bash
#!/usr/bin/env bash

set -e

conda env create -f environment.yml
conda activate biotite

mkdir -p "$CONDA_PREFIX/share/libcifpp"

curl -L --fail --retry 3 \
    -o "$CONDA_PREFIX/share/libcifpp/components.cif.gz" \
    https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz

echo
echo "Checking installation..."
python -c "import biotite; print('Biotite:', biotite.__version__)"
mkdssp --version
```

Then future installation is basically:

```bash
bash setup_biotite_env.sh
```

One caveat: `conda activate` inside a shell script can sometimes require conda initialization. The most robust version is:

```bash
conda env create -f environment.yml

conda run -n biotite bash -c '
    mkdir -p "$CONDA_PREFIX/share/libcifpp"

    curl -L --fail --retry 3 \
        -o "$CONDA_PREFIX/share/libcifpp/components.cif.gz" \
        https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz

    python -c "import biotite; print(\"Biotite:\", biotite.__version__)"
    mkdssp --version
'
```

That avoids activation entirely.

For **exact archival reproducibility on Linux**, I’d additionally save:

```bash
conda activate biotite
conda list --explicit > biotite-linux-64-spec.txt
```

That records exact builds. It can later be reconstructed with:

```bash
conda create -n biotite --file biotite-linux-64-spec.txt
```

So I’d commit these three things to your repo:

```text
environment.yml
setup_biotite_env.sh
biotite-linux-64-spec.txt
```

The first two are what humans should normally use; the explicit spec is the “recreate exactly what worked on this machine” fallback.
