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

That does not in and of itself guarantee anything about the reasonability of the backbones, so our next experiment is to calculate Phi/Psi angle distributions of our designs to assess whether we are sampling reasonable conformational space. That is the next experiment. 


# Code Review

My assessment of 4bdb42af496f7e77c92a2db85afeb77788ea204b

Broadly, I consider the code changes to be split amongst two axes: functional code changes (that which impact the execution of RFD3) and documentation changes (how to use the cylic implementation, how the args work, examples, etc). I will start by reviewing the functional code changes. 

`models/rfd3/src/rfd3/inference/input_parsing.py`

This script parses and validates inputs to build classes `DesignInputSpecification` and and the deprecated/legacy and `LegacySpecification`. It also contains various related APIs and utils. I will walkthrough the diff in stages to understand their impact. 

```sh
diff --git a/models/rfd3/src/rfd3/inference/input_parsing.py b/models/rfd3/src/rfd3/inference/input_parsing.py
index 278deb0..8dbd856 100644
--- a/models/rfd3/src/rfd3/inference/input_parsing.py
+++ b/models/rfd3/src/rfd3/inference/input_parsing.py
@@ -24,6 +24,7 @@ from pydantic import (
     BaseModel,
     ConfigDict,
     Field,
+    field_validator,
     model_validator,
 )
```
Import of the pydantic field_validator that is used later to validate i/o in `validate_cyclic_chains`.

```sh
 from rfd3.constants import (
@@ -144,6 +145,14 @@ class DesignInputSpecification(BaseModel):
     cif_parser_args: Optional[Dict[str, Any]] = Field(None, description="CIF parser arguments")
     extra: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Extra metadata to include in output (useful for logging additional info in metadata)")
     dialect: int = Field(2, description="RFdiffusion3 input dialect. 1: legacy, 2: release.")
+    cyclic_chains: list[str] | None = Field(
+        None,
+        description=(
+            "Assembled design chain IDs to use for N-to-C cyclic residue positional "
+            "encoding. Supports one complete de novo canonical peptide in dialect 2; "
+            "None or [] keeps linear encoding."
+        ),
+    )
```

Pydantic style field assignment for `cyclic_chains`. The description is accurate to the implementation, but limiting to de novo only maybe excessive. Marked for review.
 
```sh
     # ========================================================================
     # Conditioning
@@ -220,11 +229,10 @@ class DesignInputSpecification(BaseModel):
             raise FileNotFoundError(f"Output file not found at {path}")
         with open(path, "r") as f:
             data = json.load(f)
-        if "input_specification" in data:
-            spec_args = data["input_specification"]
-            return cls(**spec_args)
-        else:
+        spec_args = data.get("specification", data.get("input_specification"))
+        if spec_args is None:
             raise ValueError(f"No input specification found in json output: {path}")
+        return cls(**spec_args)
```
This is a suspicious change. I think I understand it better now. Let's bring in the whole function. 

This function is an alternative way to construct `DesignInputSpecification`.

```py
    @classmethod # decorator that let's a class use this
    def from_rfd3_out(cls, path: str):
        """Load from path to rfd3 outputs, either .cif, .cif.gz, .json or denoised / noisy trajectory files"""
        path = path.replace(".cif.gz", ".cif").replace(".cif", ".json") # cif.gz -> .cif -> .json, we get design_x.cig.gz and design_x.json in normal outputs 
        if not os.path.exists(path):
            raise FileNotFoundError(f"Output file not found at {path}")
        with open(path, "r") as f:
            data = json.load(f)
        if "input_specification" in data: # this part can't find anything (see below)
            spec_args = data["input_specification"]
            return cls(**spec_args) # in theory, constructs DesignInputSpecification
        else:
            raise ValueError(f"No input specification found in json output: {path}")
```

As an example job (run after cyclic implementation), we can look at: `experiments/007_RFD3_macrocycles/test_1/monomer/monomer_design_monomer_0_model_0.json`.

It has many json fields, but it does not have `input_specification`, it only has `specification`: 

```json
    "specification": {
        "length": "20-40",
        "extra": {
            "example": "monomer",
            "task_name": "monomer_design_monomer",
            "sampled_contig": "33P",
            "num_tokens_in": 33,
            "num_residues_in": 33,
            "num_chains": 1,
            "num_atoms": 165,
            "num_residues": 34,
            "example_id": "monomer_design_monomer_0"
        },
        "is_non_loopy": true
    },
```

It can never find that. In all of my RFD3 runs, I have never seen this error: `raise ValueError(f"No input specification found in json output: {path}")`, so I am not positive what conditions makes this loud, as I assume it should have raised this error.

So what did the agent do? 

```py
    @classmethod
    def from_rfd3_out(cls, path: str):
        """Load from path to rfd3 outputs, either .cif, .cif.gz, .json or denoised / noisy trajectory files"""
        path = path.replace(".cif.gz", ".cif").replace(".cif", ".json")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Output file not found at {path}")
        with open(path, "r") as f:
            data = json.load(f)
        spec_args = data.get("specification", data.get("input_specification")) # we prefer `specification`
        if spec_args is None:
            raise ValueError(f"No input specification found in json output: {path}")
        return cls(**spec_args)
```

The agent "fixed" this. This *may* be a good idea, but I do not know the developers intent and I have never seen this triggered. I do not think this should be part of the PR. It can be a seperate thing.

There is a chance it is important, if this is true: 

cyclic_chains saved in output JSON
        ↓
need to demonstrate it can be reloaded
        ↓
existing from_rfd3_out() cannot read current output JSON
        ↓
agent fixes from_rfd3_out()

We can review this. 
 
 ```sh
     def get_dict_to_save(self, exclude_extra: bool = False) -> dict:
         # Returns dictionary for saving (reproducible) outputs to json
@@ -371,6 +379,34 @@ class DesignInputSpecification(BaseModel):
     # Post-Validation
     # ========================================================================
 
+    @field_validator("cyclic_chains", mode="before") # This function validates the cyclic_chains field, and Pydantic should run it before Pydantic performs its normal validation/conversion of that field.
+    @classmethod # means Pydantic calls this on the class,
+    def validate_cyclic_chains(cls, value: Any) -> list[str] | None: 
+        """Validate the optional list of assembled cyclic chain IDs."""
+        if value is None:
+            return value
+        if not isinstance(value, list) or any( # is it a list? if not, raise error, "if not True or False" is what is being asked, True here is a problem
+            not isinstance(chain, str) or not chain.strip() for chain in value # are the elements str? good. if strip turns them into "" -> false, not false -> true 
+        ):
+            raise ValueError("cyclic_chains must be a list of nonempty chain IDs.")
+        if len(set(value)) != len(value):
+            raise ValueError("cyclic_chains must not contain duplicate chain IDs.")
+        if len(value) > 1:
+            raise ValueError("cyclic_chains supports at most one chain per design.")
+        return value
+
+    @model_validator(mode="after")
+    def validate_cyclic_modes(self) -> "DesignInputSpecification":
+        """Reject modes that do not support cyclic positional encoding."""
+        if self.cyclic_chains:
+            if self.dialect < 2:
+                raise ValueError("cyclic_chains requires dialect 2.")
+            if self.is_partial_diffusion:
+                raise ValueError("cyclic_chains does not support partial diffusion.")
+            if self.symmetry is not None and self.symmetry.id:
+                raise ValueError("cyclic_chains does not support symmetry.")
+        return self
+
```
I don't like posting LLM text in my personal review, but this is how I feel: 

> The cyclic field validation is internally clear and defensive. The restrictions on dialect 1 are understandable. The symmetry and partial-diffusion restrictions are left unresolved pending review of downstream assumptions; they may be technically necessary or simply conservative scope enforcement.

```sh
     @model_validator(mode="after")
     def assert_exclusivity(self):
         with validator_context("assert_exclusivity"):
@@ -515,6 +551,7 @@ class DesignInputSpecification(BaseModel):
         # Apply post-processing
         atom_array = self._append_ligand(atom_array, atom_array_input_annotated)
         atom_array = self._apply_symmetry(atom_array, atom_array_input_annotated)
+        self._validate_cyclic_chain(atom_array)
 
         # Apply globals to all tokens (including diffused)
         atom_array = self._set_origin(atom_array)
@@ -546,6 +583,33 @@ class DesignInputSpecification(BaseModel):
     # Building functions
     # ============================================================================
```

```sh 
+    def _validate_cyclic_chain(self, atom_array: AtomArray) -> None:
+        """Check only the selected chain's chemistry and de novo provenance."""
+        if not self.cyclic_chains:
+            return
+        chain_id = self.cyclic_chains[0]
+        selected = atom_array[atom_array.chain_id == chain_id]
+        if len(selected) == 0:
+            raise ValueError(f"Cyclic chain {chain_id!r} is absent from the design.") # !r just transforms A to `A`, not critical but useful
+        is_de_novo = all(
+            component.endswith("P") and component[:-1].isdigit()
+            for component in np.unique(selected.src_component)
+        )
+        annotation_names = selected.get_annotation_categories()
+        is_atomized = (
+            np.any(selected.atomize) if "atomize" in annotation_names else False
+        )
+        if (
+            not is_de_novo
+            or not np.all(np.isin(selected.res_name, STANDARD_AA))
+            or is_atomized
+            or np.any(selected.is_motif_atom_unindexed)
+        ):
+            raise ValueError(
+                f"Cyclic chain {chain_id!r} must be a complete de novo canonical "
+                "peptide, with no source-derived, atomized, or unindexed residues."
+            )
+
```
Sorry for LLM commentary again, but I do agree: 

> Understood. src_component labels generated protein blocks as e.g. 33P; this check therefore requires the entire cyclic chain to be generated de novo. I do not yet see a cyclic-RPE implementation invariant requiring this. Likely scope enforcement rather than technical necessity; marked for review.

I do not think we need the is_de_novo check. 

I also don't know if the atomization check is good, because if it is de novo, I do not believe this can happen. AtomWorks happens later. Our input here would already need to have atomize in it somehow.

> is_atomized expresses a valid downstream invariant—cyclic length assumes residue-level tokenization—but this particular validation may be redundant for the currently supported de novo canonical peptide workflow. It executes before normal pipeline atomization, and canonical amino acids are subsequently excluded from atomization unless explicitly pre-flagged. Mark for review as defensive/future-proof validation rather than obviously necessary validation.

The unindexing check *is* more reasonable, but I have not entirely drilled it down. Unindexing is user supplied and part of the design specification, there may be a way to write a CLI command that does indeed cause unindexing of the macrocycle. I am not sure. 

Conclusion so far: 

```text
is_motif_atom_unindexed
→ definitely exists before this validator
→ directly driven by a public user option

atomize
→ normally created later
→ no ordinary DesignInputSpecification field requests it
→ could already exist on a manually supplied AtomArray
```
Also, there is a warning about this sort of thing in a sense, in the docs: `atom_array_input | internal | Pre-loaded AtomArray (not recommended).`

```sh
     def _build_init(self, atom_array_input_annotated):
         # ... Fetch tokens
         indexed_tokens = (
@@ -834,6 +898,8 @@ class DesignInputSpecification(BaseModel):
     @classmethod
     def safe_init(cls, **spec_kwargs):
         if spec_kwargs.get("dialect", 2) < 2:
+            if spec_kwargs.get("cyclic_chains"):
+                raise ValueError("cyclic_chains requires dialect 2.")
             warn = (
                 "Using dialect==1, which is deprecated and will be removed in future releases. "
                 "Please update your input specification to dialect=2 and use the new schema if possible"
@@ -926,6 +992,8 @@ def create_atom_array_from_design_specification(
     **spec_kwargs,
 ) -> tuple[AtomArray, dict]:
     if int(spec_kwargs.get("dialect", 2)) < 2:
+        if spec_kwargs.get("cyclic_chains"):
+            raise ValueError("cyclic_chains requires dialect 2.")
         warn = (
             "Using dialect==1, which is deprecated and will be removed in future releases. "
             "Please update your input specification to dialect=2 and use the new schema if possible"
```

Dialect 1 does not go through the new DesignInputSpecification path where the cyclic field and its transport are implemented.

We are done with `input_parsing.py`. We've marked points for review, some of which depend on other script information we find.

`util_transforms.py`

This is a collection of RFD3-specific data-processing utilities and AtomWorks Transform implementations. Some functions are supporting helpers, while several transform classes actively construct or modify features during the inference/training pipeline.

diff --git a/models/rfd3/src/rfd3/transforms/util_transforms.py b/models/rfd3/src/rfd3/transforms/util_transforms.py
index 8b1789f..9bc739a 100644
--- a/models/rfd3/src/rfd3/transforms/util_transforms.py
+++ b/models/rfd3/src/rfd3/transforms/util_transforms.py
@@ -363,6 +363,39 @@ class RemoveTokensWithoutCorrespondingCentralAtom(Transform):
         return data
 
 ```sh
+def _resolve_cyclic_asym_ids(
+    atom_array: AtomArray,
+    token_starts: np.ndarray,
+    asym_id: np.ndarray,
+    chain_id: str,
+) -> np.ndarray:
+    """Resolve a selected peptide to its exclusive model chain identity."""
+    tokens = atom_array[token_starts]
+    selected = tokens.chain_id == chain_id
+    selected_ids = np.unique(asym_id[selected])
+    if len(selected_ids) != 1:
+        raise ValueError(
+            f"Cyclic chain {chain_id!r} must contain tokens with exactly one asym_id."
+        ) # this may be excessive in scope, but it is required in the current implementation
+    if not np.array_equal(asym_id == selected_ids[0], selected):
+        raise ValueError(
+            f"Cyclic chain {chain_id!r} must not share its asym_id with other chains."
+        )
+    selected_atoms = atom_array[atom_array.chain_id == chain_id]
+    selected_residue_indices = tokens.within_chain_res_idx[selected]
+    if (
+        np.any(selected_atoms.atomize)
+        or not np.all(np.isin(selected_atoms.res_name, STANDARD_AA))
+        or np.any(selected_atoms.is_motif_atom_unindexed)
+        or not np.all(np.diff(selected_residue_indices) == 1) # is it contigious 
+    ):
+        raise ValueError(
+            f"Cyclic chain {chain_id!r} requires one canonical, nonatomized, indexed "
+            "token per consecutive residue."
+        )
+    return selected_ids.astype(np.int64)
+
+
```
The user-selected chain must correspond to exactly one model chain identity.

The implementation currently supports exactly one cyclic chain and, more deeply, requires that the selected public chain ID resolve to exactly one model asym_id. This makes the implementation robust for the intended monomer/single-binder workflows but intentionally excludes otherwise scientifically plausible multi-macrocycle cases. It is not yet clear whether a public chain ID that maps to multiple chain instances should be considered ambiguous or should naturally mean “apply cyclic encoding to all instances.” This is a repository/API design decision worth asking maintainers rather than guessing.

Also, this maybe too defensive but its impossible to know without an understanding of what happens in atomworks. At the point EncodeAF3TokenLevelFeatures runs, does every chain_id necessarily correspond to exactly one asym_id for every input RFD3 permits? 

This is a defensive one-to-one mapping check. It may be redundant for ordinary monomer/single-binder workflows, but AtomWorks does not generally guarantee that a public chain_id identifies exactly one chain instance. Since failure of this assumption could silently select the wrong RPE target, keeping the check is reasonable.

```sh
 class EncodeAF3TokenLevelFeatures(Transform):
     def __init__(
         self,
@@ -409,6 +442,13 @@ class EncodeAF3TokenLevelFeatures(Transform):
         # ... (within chain entity)
         sym_name, sym_id = get_within_entity_idx(token_level_array, level="pn_unit")
 
+        # Validate residue identity before sequence masking replaces names with GAP.
+        cyclic_chains = data.get("specification", {}).get("cyclic_chains") or []
+        if cyclic_chains:
+            data.setdefault("feats", {})["cyclic_asym_ids"] = _resolve_cyclic_asym_ids(
+                atom_array, token_starts, asym_id, cyclic_chains[0]
+            )
+
         # ... molecule type
         _aa_like_res_names = self.sequence_encoding.all_res_names[
             self.sequence_encoding.is_aa_like
```

`models/rfd3/src/rfd3/model/layers/blocks.py`

diff --git a/models/rfd3/src/rfd3/model/layers/blocks.py b/models/rfd3/src/rfd3/model/layers/blocks.py
index 9963d6e..da415c6 100644
--- a/models/rfd3/src/rfd3/model/layers/blocks.py
+++ b/models/rfd3/src/rfd3/model/layers/blocks.py
@@ -274,6 +274,32 @@ class SimpleRecycler(nn.Module):
         return S_I, Z_II
 
 ```sh
+def _cyclic_residue_offsets(
+    offsets: torch.Tensor, asym_id: torch.Tensor, cyclic_asym_ids: torch.Tensor
+) -> torch.Tensor:
+    """Wrap intrachain offsets for one canonical peptide, preserving half-ring ties.
+
+    Args:
+        offsets: Signed residue differences, with shape [I, I].
+        asym_id: Chain identities, with shape [I].
+        cyclic_asym_ids: The selected chain identity, with shape [1]. Input
+            validation guarantees consecutive residues and one token per residue.
+
+    Returns:
+        Residue offsets with only the selected intrachain pairs wrapped.
+    """
+    selected = asym_id == cyclic_asym_ids[0]
+    length = selected.sum()
+    wrapped = torch.where(
+        2 * offsets > length,
+        offsets - length,
+        torch.where(2 * offsets < -length, offsets + length, offsets),
+    )
+    return torch.where(
+        selected.unsqueeze(-1) & selected.unsqueeze(-2), wrapped, offsets
+    )
+
+
```


```sh
 class RelativePositionEncodingWithIndexRemoval(nn.Module):
     """
     Usual RPE but utilizes `is_motif_atom_3d_unindexed` to ensure within-chain position is spoofed.
@@ -295,12 +321,18 @@ class RelativePositionEncodingWithIndexRemoval(nn.Module):
     def forward(self, f):
         b_samechain_II = f["asym_id"].unsqueeze(-1) == f["asym_id"].unsqueeze(-2)
         b_same_entity_II = f["entity_id"].unsqueeze(-1) == f["entity_id"].unsqueeze(-2)
+        residue_offsets = f["residue_index"].unsqueeze(-1) - f[
+            "residue_index"
+        ].unsqueeze(-2)
+        cyclic_asym_ids = f.get("cyclic_asym_ids")
+        if cyclic_asym_ids is not None and cyclic_asym_ids.numel() > 0:
+            residue_offsets = _cyclic_residue_offsets(
+                residue_offsets, f["asym_id"], cyclic_asym_ids
+            )
         d_residue_II = torch.where(
             b_samechain_II,
             torch.clip(
-                f["residue_index"].unsqueeze(-1)
-                - f["residue_index"].unsqueeze(-2)
-                + self.r_max,
+                residue_offsets + self.r_max,
                 0,
                 2 * self.r_max,
             ),
```

This is a REAL divergence. My proposal here was a simple if/else. if cyclic is requested, do cyclic stuff, otherwise do the normal linear encoding verbatim. This is how RF3 handles it. 

Note:
b_... → boolean. So b_samechain_II is a Boolean mask.
d_... → usually difference/distance/offset. For example d_residue_II is the pairwise residue-index difference after processing.
_I → one token/residue dimension, conceptually shape [N].
_II → two token/residue dimensions, conceptually shape [N, N]: an i,j pair representation.

torch.where(condition, x, y)

IF 2*offset > length:
    use offset - length
ELSE:
    evaluate the inner torch.where

if 2 * offset > length:
    wrapped = offset - length
elif 2 * offset < -length:
    wrapped = offset + length
else:
    wrapped = offset


The first major deviation is the introduction of a `residue_offset` term. These values were present prior to the implementation, but not assigned to a variable. 

It is defined: 

```sh
        residue_offsets = f["residue_index"].unsqueeze(-1) - f[
            "residue_index"
        ].unsqueeze(-2)
```

This makes a row value - column value 2D table. Note that here, it is not confined as an intra-chain calculation. Otherwise, it behaves the same. Here is an example table for a monomer design of 6 canonical AAs, followed by an analogous table for a design of 5 canonical AAs.

    i_1	i_2	I_3	i_4	i_5	i_6
j_1	0	-1	-2	-3	-4	-5
j_2	1	0	-1	-2	-3	-4
j_3	2	1	0	-1	-2	-3
j_4	3	2	1	0	-1	-2
j_5	4	3	2	1	0	-1
j_6	5	4	3	2	1	0

    i_1	i_2	I_3	i_4	i_5
j_1	0	-1	-2	-3	-4
j_2	1	0	-1	-2	-3
j_3	2	1	0	-1	-2
j_4	3	2	1	0	-1
j_5	4	3	2	1	0

The above (and variations due to the amount of chains specified) will run regardless of our divergence. 

```sh
cyclic_asym_ids = f.get("cyclic_asym_ids")
        if cyclic_asym_ids is not None and cyclic_asym_ids.numel() > 0:
```

We check that we have a cyclic_asm_ids assignment. It can not be `None` and `numel()` (the amount of elements in the tensor) must be > 0. 

If that passes, we enter the cyclic encoding function to update the residue offsets shown in the table. 

```sh
        cyclic_asym_ids = f.get("cyclic_asym_ids")
        if cyclic_asym_ids is not None and cyclic_asym_ids.numel() > 0:
            residue_offsets = _cyclic_residue_offsets(
                residue_offsets, f["asym_id"], cyclic_asym_ids
            )
```
We enter `_cyclic_residue_offsets()`. 

```sh
def _cyclic_residue_offsets(
    offsets: torch.Tensor, asym_id: torch.Tensor, cyclic_asym_ids: torch.Tensor
) -> torch.Tensor:
    """Wrap intrachain offsets for one canonical peptide, preserving half-ring ties.

    Args:
        offsets: Signed residue differences, with shape [I, I].
        asym_id: Chain identities, with shape [I].
        cyclic_asym_ids: The selected chain identity, with shape [1]. Input
            validation guarantees consecutive residues and one token per residue.

    Returns:
        Residue offsets with only the selected intrachain pairs wrapped.
    """
    selected = asym_id == cyclic_asym_ids[0] 
    length = selected.sum()
    wrapped = torch.where(
        2 * offsets > length,
        offsets - length,
        torch.where(2 * offsets < -length, offsets + length, offsets),
    )
    return torch.where(
        selected.unsqueeze(-1) & selected.unsqueeze(-2), wrapped, offsets
    )
```
We create a boolean mask from the equality check of tensor asym_id to tensor cyclic_asym_ids[0], True corresponds to residue/token of the macrocycle. Our prior parsing guarantees consecutive residues and one token per residue. We get the length of our desired macrocycle by summing True states in the boolean mask in `selected`. 

We begin the wrapping process. The nested torch.where() can hard to read. Remember: `torch.where(CONDITION, x, y)`. We assess if `2 * offsets > length` is True for an element in the tensor, if so, we replace that element with `offsets - length`. If `2 * offsets > length` is False, we enter the second torch.where operation: `torch.where(2 * offsets < -length, offsets + length, offsets)`. Here, we assess the condition to find the value we will assign to `y` of the first `torch.where()`. If `2 * offsets < -length` is True, we will return `offsets + length`, if False: `offsets`. 

Let's work through an example with a 6-member chain.

Offsets is: 

    i_1	i_2	I_3	i_4	i_5	i_6
j_1	0	-1	-2	-3	-4	-5
j_2	1	0	-1	-2	-3	-4
j_3	2	1	0	-1	-2	-3
j_4	3	2	1	0	-1	-2
j_5	4	3	2	1	0	-1
j_6	5	4	3	2	1	0

The diagonal will remain `0` throughout. i_1/j_2 is scalar `1`. We evaluate the CONDITION, `2 * offsets > length`. We check if (1 * 2 > 6). It is not, so we enter the torch.where(). We evaluate the condition `torch.where(2 * offsets < -length`. We check if (2 * 1 < -6). It is not. So we return offsets as is, `1`. The same process will occur for i_1/j_3, and we return 2. At i_1/j_4, we are at the halfway point, we still do not satisfy either condition so we return `3`. At i_1/j_5, we do satisfy (4 * 2 > 6). We take the shortest way around the ring, 4-6 = -2, we return -2. At i_1/j_6, we again satisfy that condition: ( 2 * 5 > 6). We take the shortest route, 5-6 = -1.

Let's move to i_2. At i_2/j_1, we have `-1` in the `offset`. Evaluate: (2 * (-1) > 6) is False. Evaluate (2 * (-1) < -6) is False. Return `-1`. Let us continue moving across the i direction. At i_3/j_1: Evaluate: (2 * (-2) > 6) is False. Evaluate (2 * (-2) < -6) is False. Return `-2`. Same at i_4/j_1, return `-3`. At i_5/j_1: Evaluate (2 * (-4) > 6), False -> (2 * -4 < -6), True. So we return `offsets + length`: -4 + 6 = `2`, return 2. This will happen again at i_6/j_1, we return `1`. 

We can go through every i and j position like this. We get: 

    i_1	i_2	I_3	i_4	i_5	i_6
j_1	0	-1	-2	-3	2	1
j_2	1	0	-1	-2	-3	2
j_3	2	1	0	-1	-2	-3
j_4	3	2	1	0	-1	-2
j_5	-2	3	2	1	0	-1
j_6	-1	-2	3	2	1	0


We now have `wrapped`. We will only wrap the conditions where `selected.unsqueeze(-1) & selected.unsqueeze(-2)` is true, recall that: `selected = asym_id == cyclic_asym_ids[0]`. 


**I went through the rest of this elsewhere**

Tests: 

`models/rfd3/tests/test_rfd3_cyclic_encoding.py`

- I removed the MPS test, but otherwise this code is broadly fine. I think it relies on magic numbers too strongly and can be difficult to read, but the implementation of the tests is correct. 

`models/rfd3/tests/test_rfd3_cyclic_inputs.py`

all cyclic tests that actually execute to completion pass; the two expected failures are hitting an existing Foundry inference/annotation problem that also occurs with cyclic encoding disabled.

```sh
def test_one_residue_and_sampled_length_are_supported():
    one_residue = DesignInputSpecification(length="1", cyclic_chains=["A"])
    assert len(get_token_starts(one_residue.build())) == 1

    sampled = DesignInputSpecification(length="2-4", cyclic_chains=["A"])
    sampled_array = sampled.build()
    assert 2 <= len(get_token_starts(sampled_array)) <= 4
```

one residue arguably *shouldn't* be supported, so this test is backwards. sampling still being possible is more of an upstream problem. 

```sh
@pytest.mark.parametrize("cyclic_chains", [None, ["A"]])
def test_fixed_cofactor_survives_pipeline_with_or_without_cyclic_encoding(
    cyclic_chains, inference_pipeline
):
    spec = DesignInputSpecification(
        input=str(COFACTOR_INPUT),
        contig="3,/0,B1-5,C1-1",
        select_fixed_atoms=True,
        cyclic_chains=cyclic_chains,
    )
    pipeline_input = spec.to_pipeline_input("cofactor")
    assert "NAG" in pipeline_input["atom_array"].res_name
    try:
        output = inference_pipeline(pipeline_input)
    except AttributeError as error:
        if "is_motif_atom" not in str(error):
            raise
        pytest.xfail(
            "Pre-existing ligand reference annotation failure, reproduced without "
            "cyclic encoding"
        )
    assert "NAG" in output["atom_array"].res_name
    if cyclic_chains:
        _assert_cyclic_rpe(output, "A")
```

This is a question of if they want it or not, I don't think they do, but this one is easy to read and not harmful. 

```sh
def test_unsupported_modes_are_rejected_at_public_dispatches():
    with pytest.raises(ValueError, match="dialect 2"):
        DesignInputSpecification.safe_init(dialect=1, length="3", cyclic_chains=["A"])
    with pytest.raises(ValueError, match="dialect 2"):
        create_atom_array_from_design_specification(
            dialect=1, length="3", cyclic_chains=["A"]
        )
    with pytest.raises(ValidationError, match="partial diffusion"):
        DesignInputSpecification(input=str(TARGET), partial_t=1.0, cyclic_chains=["A"])
    with pytest.raises(ValidationError, match="symmetry"):
        DesignInputSpecification(
            length="3",
            symmetry={"id": "C2", "is_symmetric_motif": False},
            cyclic_chains=["A"],
        )
```

Same. 

```sh
def test_build_rejects_absent_source_derived_and_mixed_chains():
    with pytest.raises(ValueError, match="absent"):
        DesignInputSpecification(length="3", cyclic_chains=["B"]).build()

    source = DesignInputSpecification(
        input=str(TARGET),
        contig="E6-10",
        select_fixed_atoms=False,
        select_unfixed_sequence=True,
        cyclic_chains=["A"],
    )
    with pytest.raises(ValueError, match="complete de novo"):
        source.build()

    mixed = DesignInputSpecification(
        input=str(TARGET), contig="3,E6-10", cyclic_chains=["A"]
    )
    with pytest.raises(ValueError, match="complete de novo"):
        mixed.build()

    atomized = DesignInputSpecification(length="3", cyclic_chains=["A"])
    atomized_array = atomized.build()
    atomized_array.set_annotation(
        "atomize", [True] + [False] * (len(atomized_array) - 1)
    )
    with pytest.raises(ValueError, match="atomized"):
        atomized._validate_cyclic_chain(atomized_array)
```
["B"] with length="3"    → requested chain doesn't exist
"E6-10"                  → cyclic chain exists, but is source-derived
"3,E6-10"                → cyclic chain exists, but is mixed


“For now, cyclic chains must consist entirely of normal residue-level tokens. If any part of the requested cyclic chain is atomized, reject it rather than pretending the RPE semantics are well defined.”

```sh
@pytest.mark.parametrize("metadata_key", ["specification", "input_specification"])
def test_output_metadata_round_trips_cyclic_request(tmp_path, metadata_key):
    spec = DesignInputSpecification(length="3", cyclic_chains=["A"])
    output = tmp_path / "design.json"
    output.write_text(json.dumps({metadata_key: spec.get_dict_to_save()}))

    restored = DesignInputSpecification.from_rfd3_out(str(output))
    assert restored.cyclic_chains == ["A"]
```

Again, did `DesignInputSpecification.from_rfd3_out()` ever work? I don't think so, so no need to continute supporting two keys. It depends on how they view my change there, because I made it work. 

