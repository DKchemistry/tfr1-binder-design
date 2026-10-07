reinitialize 

python
from pathlib import Path
import re
from pymol import cmd

OBJ = "oracle"

files = list(Path("oracle").glob("rfdiffusion_*/round_1/prediction.pdb"))
files.sort(
    key=lambda p: int(re.search(r"rfdiffusion_(\d+)", str(p)).group(1))
)

if not files:
    raise RuntimeError(
        "No files found matching oracle/rfdiffusion_*/round_1/prediction.pdb"
    )

cmd.delete(OBJ)
cmd.set("all_states", 0)

for state, pdb in enumerate(files, start=1):
    match = re.search(r"rfdiffusion_(\d+)", str(pdb))
    design_id = int(match.group(1))

    oracle_name = f"oracle_{design_id}"

    print(f"PyMOL state {state}: {oracle_name} <- {pdb}")

    cmd.load(
        str(pdb),
        OBJ,
        state=state,
        discrete=1,
        zoom=0,
    )

    # State title matches the actual source directory number exactly
    cmd.set_title(
        OBJ,
        state,
        oracle_name,
    )

# Recalculate secondary structure independently for every state
for state in range(1, len(files) + 1):
    cmd.dss(
        OBJ,
        state=state,
    )

cmd.rebuild()

# Start on oracle_0
cmd.set("state", 1)

# Fit all structures into view
cmd.zoom(OBJ)

print()
print(f"Loaded {len(files)} Oracle structures into object '{OBJ}'")
print("State mapping:")
for state, pdb in enumerate(files, start=1):
    design_id = int(re.search(r"rfdiffusion_(\d+)", str(pdb)).group(1))
    print(f"  state {state:3d} -> oracle_{design_id}")

python end

hide everything, oracle
show cartoon, oracle
cartoon automatic, oracle