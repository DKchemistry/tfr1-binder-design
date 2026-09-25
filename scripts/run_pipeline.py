#!/usr/bin/env python3

import argparse
import os
import shlex
import subprocess
import tomllib
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent


def resolve_path(value):
    path = Path(value).expanduser()

    if path.is_absolute():
        return path.resolve()

    return (REPO_ROOT / path).resolve()


def run_command(command, cwd=None, env=None):
    command = [str(item) for item in command]

    print(f"\n$ {shlex.join(command)}")

    subprocess.run(
        command,
        check=True,
        cwd=cwd,
        env=env,
    )


def load_config(config_path):
    with open(config_path, "rb") as handle:
        return tomllib.load(handle)


def run_tensors(config, paths):
    print("\n" + "=" * 70)
    print("1. Interface tensors")
    print("=" * 70)

    tensor_config = config["tensors"]

    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["tensors"],
        "python",
        str(REPO_ROOT / "scripts/interface_tensors/make_interface_tensor.py"),

        "--input_pdb",
        str(paths["target_pdb"]),

        "--out_dir",
        str(paths["tensor_dir"]),

        "--binderlen",
        str(tensor_config["binder_length"]),

        "--target_adj",
        tensor_config["target_adj"],

        "--binder_ss",
        tensor_config["binder_ss"],

        "--binder_ss_len",
        str(tensor_config["binder_ss_len"]),
    ]

    run_command(command)


def run_rfdiffusion(config, paths):
    print("\n" + "=" * 70)
    print("2. RFdiffusion")
    print("=" * 70)

    rfd_config = config["rfdiffusion"]

    hotspots = ",".join(rfd_config["hotspots"])

    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["rfdiffusion"],
        "python",
        "scripts/run_inference.py",

        "--config-name",
        "base",

        f"inference.output_prefix={paths['rfdiffusion_dir'] / 'rfdiffusion'}",

        f"inference.num_designs={rfd_config['num_designs']}",

        f"inference.input_pdb={paths['target_pdb']}",

        f"contigmap.contigs=[{rfd_config['contig']}]",

        "inference.model_runner=ScaffoldedSampler",

        "scaffoldguided.scaffoldguided=True",

        f"scaffoldguided.scaffold_dir={paths['tensor_dir']}",

        f"scaffoldguided.scaffold_list={paths['tensor_list']}",

        f"inference.cyclic={rfd_config['cyclic']}",

        f"inference.cyc_chains={rfd_config['cyc_chains']}",

        f"diffuser.T={rfd_config['diffusion_steps']}",

        f"ppi.hotspot_res=[{hotspots}]",
    ]

    run_command(
        command,
        cwd=paths["rfdiffusion_root"],
    )


def run_distal_site(config, paths):
    print("\n" + "=" * 70)
    print("3. Distal-site selection")
    print("=" * 70)

    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["biopython"],
        "python",
        str(REPO_ROOT / "scripts/select_distal_site.py"),

        str(paths["rfdiffusion_dir"]),

        "--peptide-chain",
        config["distal_site"]["peptide_chain"],

        "--output",
        str(paths["distal_csv"]),
    ]

    run_command(command)


def run_mpnn_relax(config, paths):
    print("\n" + "=" * 70)
    print("4. ProteinMPNN / RosettaRelax")
    print("=" * 70)

    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["biopython"],
        "python",
        str(REPO_ROOT / "scripts/run_mpnn_relax_all.py"),

        "--pdb-dir",
        str(paths["rfdiffusion_dir"]),

        "--scores",
        str(paths["distal_csv"]),

        "--output-dir",
        str(paths["mpnn_relax_dir"]),

        "--xml",
        str(paths["relax_xml"]),

        "--proteinmpnn-dir",
        str(paths["proteinmpnn_root"]),

        "--rounds",
        str(config["mpnn_relax"]["rounds"]),
    ]

    run_command(command)


def run_oracle(config, paths):
    print("\n" + "=" * 70)
    print("5. AfCyc oracle")
    print("=" * 70)

    oracle_config = config["oracle"]

    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["biopython"],
        "python",
        str(REPO_ROOT / "scripts/run_afcyc_all.py"),

        "--design-dir",
        str(paths["mpnn_relax_dir"]),

        "--output-dir",
        str(paths["oracle_dir"]),

        "--target-pdb",
        str(paths["target_pdb"]),

        "--target-chain",
        oracle_config["target_chain"],

        "--params",
        str(paths["alphafold_params"]),

        "--afcyc-env",
        config["envs"]["afcyc"],

        "--recycles",
        str(oracle_config["recycles"]),

        "--seed",
        str(oracle_config["seed"]),
    ]

    env = os.environ.copy()

    # Necessary on memory-constrained GPUs such as the current 6 GB test GPU.
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

    run_command(
        command,
        env=env,
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Pipeline TOML configuration file.",
    )

    args = parser.parse_args()

    config_path = args.config.expanduser().resolve()
    config = load_config(config_path)

    project_dir = resolve_path(config["project_dir"])
    target_pdb = resolve_path(config["target_pdb"])

    tensor_dir = (
        project_dir
        / "tensors"
        / config["tensors"]["output_name"]
    )

    paths = {
        "project_dir": project_dir,
        "target_pdb": target_pdb,

        "tensor_dir": tensor_dir,
        "tensor_list": tensor_dir.with_suffix(".txt"),

        "rfdiffusion_dir": project_dir / "rfdiffusion",
        "distal_csv": project_dir / "distal_site" / "distal_site.csv",
        "mpnn_relax_dir": project_dir / "mpnn_relax",
        "oracle_dir": project_dir / "oracle",

        "rfdiffusion_root": resolve_path(
            config["paths"]["rfdiffusion"]
        ),

        "proteinmpnn_root": resolve_path(
            config["paths"]["proteinmpnn"]
        ),

        "alphafold_params": resolve_path(
            config["paths"]["alphafold_params"]
        ),

        "relax_xml": resolve_path(
            config["paths"]["relax_xml"]
        ),
    }

    if not target_pdb.is_file():
        raise FileNotFoundError(
            f"Target PDB not found: {target_pdb}"
        )

    for directory in [
        project_dir / "tensors",
        paths["rfdiffusion_dir"],
        paths["distal_csv"].parent,
        paths["mpnn_relax_dir"],
        paths["oracle_dir"],
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    print("\n" + "=" * 70)
    print("Pipeline configuration")
    print("=" * 70)
    print(f"Config:      {config_path}")
    print(f"Project:     {project_dir}")
    print(f"Target PDB:  {target_pdb}")

    run_tensors(config, paths)
    run_rfdiffusion(config, paths)
    run_distal_site(config, paths)
    run_mpnn_relax(config, paths)
    run_oracle(config, paths)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Project: {project_dir}")


if __name__ == "__main__":
    main()