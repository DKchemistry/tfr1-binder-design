import gzip
import tempfile
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx


SEARCH_DIR = Path("../test_2/macrocycle_monomer_10K")
MIN_HELIX_LENGTH = 4


def longest_run(sse, target="a"):
    best = current = 0

    for x in sse:
        if x == target:
            current += 1
            best = max(best, current)
        else:
            current = 0

    return best


for path in sorted(SEARCH_DIR.glob("*.cif.gz")):
    try:
        with gzip.open(path, "rt") as f:
            cif = pdbx.CIFFile.read(f)

        atoms = pdbx.get_structure(cif, model=1)

        sse = struc.annotate_sse(atoms)
        helix_len = longest_run(sse, "a")

        if helix_len >= MIN_HELIX_LENGTH:
            print(f"FOUND: {path}")
            print(f"alpha helix: {helix_len} residues")
            break

    except Exception as e:
        print(f"ERROR {path}: {e}")

else:
    print("No alpha helix found.")