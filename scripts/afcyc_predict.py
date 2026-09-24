#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import numpy as np
import jax
import jax.numpy as jnp

# Some ColabDesign code paths still reference jax.tree_map, which was moved
# to jax.tree_util.tree_map in newer JAX releases.
if not hasattr(jax, "tree_map"):
    jax.tree_map = jax.tree_util.tree_map

# jnp.clip(a_min=..., a_max=...) -> jnp.clip(min=..., max=...)
#
# Older AlphaFold/ColabDesign code uses NumPy's historical a_min/a_max
# keyword names. Current JAX uses min/max instead.
_jax_clip = jnp.clip


def _clip_compat(
    arr=None,
    a_min=None,
    a_max=None,
    *,
    min=None,
    max=None,
):
    if a_min is not None and min is not None:
        raise TypeError("Specify only one of a_min or min")

    if a_max is not None and max is not None:
        raise TypeError("Specify only one of a_max or max")

    lower = min if min is not None else a_min
    upper = max if max is not None else a_max

    return _jax_clip(
        arr,
        min=lower,
        max=upper,
    )


jnp.clip = _clip_compat


# Import ColabDesign only after installing the compatibility shims.
from colabdesign import mk_afdesign_model


def add_cyclic_binder_offset(model, offset_type=2):
    """
    Apply AfCycDesign-style cyclic relative positional offsets to the binder.

    For protocol="binder":
      - target-target offsets remain unchanged
      - target-binder offsets remain unchanged
      - binder-binder offsets become cyclic

    offset_type=2 matches the AfCycDesign notebook implementation.
    """

    def cyclic_offset(length):
        i = np.arange(length)

        # Represent each residue both in the original ring and one ring-length
        # away, then find the shortest sequence-space distance.
        ij = np.stack([i, i + length], axis=-1)

        offset = i[:, None] - i[None, :]

        c_offset = np.abs(
            ij[:, None, :, None] - ij[None, :, None, :]
        ).min(axis=(2, 3))

        if offset_type == 1:
            pass

        elif offset_type >= 2:
            wrapped_is_shorter = c_offset < np.abs(offset)
            c_offset[wrapped_is_shorter] = -c_offset[wrapped_is_shorter]

        if offset_type == 3:
            idx = np.abs(c_offset) > 2
            c_offset[idx] = (
                32 * c_offset[idx] / np.abs(c_offset[idx])
            )

        return c_offset * np.sign(offset)

    if model.protocol != "binder":
        raise ValueError(
            f"Expected protocol='binder', got protocol={model.protocol!r}"
        )

    residue_index = np.asarray(model._inputs["residue_index"])

    offset = np.array(
        residue_index[:, None] - residue_index[None, :]
    )

    c_offset = cyclic_offset(model._binder_len)

    target_len = model._target_len

    offset[
        target_len:,
        target_len:
    ] = c_offset

    model._inputs["offset"] = offset


def validate_sequence(sequence):
    allowed = set("ACDEFGHIKLMNPQRSTVWY")

    sequence = sequence.strip().upper()

    invalid = sorted(set(sequence) - allowed)

    if invalid:
        raise ValueError(
            "Sequence contains unsupported residue(s): "
            + ", ".join(invalid)
        )

    if not sequence:
        raise ValueError("Sequence is empty.")

    return sequence


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Predict a target + cyclic peptide complex using "
            "ColabDesign/AfCycDesign-style cyclic positional offsets."
        )
    )

    parser.add_argument(
        "--pdb",
        required=True,
        help="Input PDB containing the target structure.",
    )

    parser.add_argument(
        "--target_chain",
        default="A",
        help="Target chain in the input PDB. Default: A",
    )

    parser.add_argument(
        "--sequence",
        required=True,
        help="Cyclic binder amino-acid sequence.",
    )

    parser.add_argument(
        "--params",
        default="~/alphafold",
        help=(
            "AlphaFold data directory containing params/. "
            "Default: ~/alphafold"
        ),
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output directory.",
    )

    parser.add_argument(
        "--recycles",
        type=int,
        default=1,
        help="Number of AlphaFold recycles. Default: 1",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed. Default: 0",
    )

    args = parser.parse_args()

    # -----------------------------------------------------------------------
    # Resolve inputs
    # -----------------------------------------------------------------------

    pdb = Path(args.pdb).expanduser().resolve()
    params = Path(args.params).expanduser().resolve()
    out = Path(args.out).expanduser().resolve()

    if not pdb.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb}")

    if not params.is_dir():
        raise FileNotFoundError(
            f"AlphaFold data directory not found: {params}"
        )

    if not (params / "params").is_dir():
        raise FileNotFoundError(
            f"Expected AlphaFold parameter directory at: "
            f"{params / 'params'}"
        )

    if args.recycles < 0:
        raise ValueError("--recycles must be >= 0")

    sequence = validate_sequence(args.sequence)

    out.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Report configuration
    # -----------------------------------------------------------------------

    print("AfCyc local prediction")
    print("----------------------")
    print(f"JAX devices:      {jax.devices()}")
    print(f"Target PDB:       {pdb}")
    print(f"Target chain:     {args.target_chain}")
    print(f"Binder sequence:  {sequence}")
    print(f"Binder length:    {len(sequence)}")
    print(f"AlphaFold params: {params}")
    print(f"Model:            model_1_ptm")
    print(f"Recycles:         {args.recycles}")
    print(f"Seed:             {args.seed}")
    print(f"Output directory: {out}")
    print()

    # -----------------------------------------------------------------------
    # Initialize AlphaFold/ColabDesign
    # -----------------------------------------------------------------------
    #
    # Deliberately load only model_1_ptm for this proof-of-concept run.
    # This avoids loading all available AF2 parameter sets into memory.
    #
    model = mk_afdesign_model(
        protocol="binder",
        data_dir=str(params),
        model_names=["model_1_ptm"],
    )

    # -----------------------------------------------------------------------
    # Prepare target + binder
    # -----------------------------------------------------------------------
    #
    # Only the target chain is read from the supplied PDB.
    #
    # binder_len tells ColabDesign to append a new binder sequence of that
    # length. The RFdiffusion coordinates for chain B therefore are NOT
    # supplied to AlphaFold as a binder template.
    #
    model.prep_inputs(
        pdb_filename=str(pdb),
        target_chain=args.target_chain,
        binder_len=len(sequence),
    )

    print(
        f"Prepared target length: {model._target_len} residues"
    )
    print(
        f"Prepared binder length: {model._binder_len} residues"
    )

    if model._binder_len != len(sequence):
        raise RuntimeError(
            "Prepared binder length does not match supplied sequence."
        )

    # -----------------------------------------------------------------------
    # Apply AfCyc cyclic positional encoding
    # -----------------------------------------------------------------------

    add_cyclic_binder_offset(
        model,
        offset_type=2,
    )

    print("Applied cyclic binder positional offsets.")
    print("Starting AlphaFold prediction...")
    print()

    # -----------------------------------------------------------------------
    # Predict
    # -----------------------------------------------------------------------

    aux = model.predict(
        seq=sequence,
        models=["model_1_ptm"],
        num_recycles=args.recycles,
        sample_models=False,
        dropout=False,
        seed=args.seed,
        return_aux=True,
        verbose=True,
    )

    # -----------------------------------------------------------------------
    # Output paths
    # -----------------------------------------------------------------------

    pdb_out = out / "prediction.pdb"
    pae_out = out / "pae.npy"
    offset_out = out / "cyclic_offset.npy"
    metrics_out = out / "metrics.json"

    # -----------------------------------------------------------------------
    # Save predicted structure
    # -----------------------------------------------------------------------

    model.save_current_pdb(str(pdb_out))

    # Raw AlphaFold PAE matrix, in Angstroms.
    pae = np.asarray(aux["pae"])
    np.save(pae_out, pae)

    # Save the actual relative-position matrix used so we can verify the
    # cyclic encoding later if necessary.
    np.save(
        offset_out,
        np.asarray(model._inputs["offset"]),
    )

    # -----------------------------------------------------------------------
    # Metrics
    # -----------------------------------------------------------------------
    #
    # In binder mode ColabDesign's logged pLDDT is calculated over the binder.
    #
    # i_pae is ColabDesign's normalized interface PAE:
    #
    #     mean interface PAE / 31 Å
    #
    # so lower is better.
    #
    log = aux["log"]

    plddt_binder = float(log["plddt"])
    i_pae_normalized = float(log["i_pae"])
    ptm = float(log["ptm"])
    i_ptm = float(log["i_ptm"])

    metrics = {
        "target_pdb": str(pdb),
        "target_chain": args.target_chain,
        "target_length": int(model._target_len),
        "sequence": sequence,
        "binder_length": len(sequence),
        "cyclic_offset_type": 2,
        "model": "model_1_ptm",
        "recycles": args.recycles,
        "seed": args.seed,
        "plddt_binder": plddt_binder,
        "plddt_binder_100": plddt_binder * 100.0,
        "i_pae_normalized": i_pae_normalized,
        "i_pae_angstrom_equivalent": i_pae_normalized * 31.0,
        "ptm": ptm,
        "i_ptm": i_ptm,
    }

    with open(metrics_out, "w") as handle:
        json.dump(metrics, handle, indent=2)

    print()
    print("Prediction complete.")
    print("--------------------")
    print(f"Binder pLDDT:       {plddt_binder:.4f}")
    print(
        f"Binder pLDDT ×100:  "
        f"{plddt_binder * 100.0:.2f}"
    )
    print(
        f"Normalized iPAE:    "
        f"{i_pae_normalized:.4f}"
    )
    print(
        f"iPAE ×31 Å:         "
        f"{i_pae_normalized * 31.0:.2f} Å"
    )
    print(f"pTM:                {ptm:.4f}")
    print(f"ipTM:               {i_ptm:.4f}")

    print()
    print("Saved:")
    print(f"  {pdb_out}")
    print(f"  {pae_out}")
    print(f"  {offset_out}")
    print(f"  {metrics_out}")


if __name__ == "__main__":
    main()