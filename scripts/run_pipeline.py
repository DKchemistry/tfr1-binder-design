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
    subprocess.run(command, check=True, cwd=cwd, env=env)


def load_config(config_path):
    with open(config_path, "rb") as handle:
        return tomllib.load(handle)


def print_stage(number, title):
    print("\n" + "=" * 70)
    print(f"{number}. {title}")
    print("=" * 70)


def run_tensors(config, paths):
    print_stage(1, "Interface tensors")
    tensor_config = config["tensors"]
    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["tensors"],
        "python", str(SCRIPT_DIR / "interface_tensors/make_interface_tensor.py"),
        "--input_pdb", str(paths["target_pdb"]),
        "--out_dir", str(paths["tensor_dir"]),
        "--binderlen", str(tensor_config["binder_length"]),
        "--target_adj", tensor_config["target_adj"],
        "--binder_ss", tensor_config["binder_ss"],
        "--binder_ss_len", str(tensor_config["binder_ss_len"]),
    ]
    run_command(command)


def run_rfdiffusion(config, paths):
    print_stage(2, "RFdiffusion")
    rfd_config = config["rfdiffusion"]
    hotspots = ",".join(rfd_config["hotspots"])
    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["rfdiffusion"],
        "python", "scripts/run_inference.py",
        "--config-name", "base",
        f"inference.output_prefix={paths['generated_backbone_dir'] / 'rfdiffusion'}",
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
    run_command(command, cwd=paths["rfdiffusion_root"])


def run_distal_site(config, paths):
    print_stage(3, "Distal-site selection")
    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["biopython"],
        "python", str(SCRIPT_DIR / "select_distal_site.py"),
        str(paths["backbone_dir"]),
        "--peptide-chain", config["distal_site"]["peptide_chain"],
        "--output", str(paths["generated_distal_csv"]),
    ]
    run_command(command)


def run_inverse_design(config, paths):
    design_config = config["mpnn_relax"]
    method = design_config.get("method", "proteinmpnn")
    relax = design_config.get("relax", True)
    design_chain = design_config.get(
        "design_chain",
        config.get("distal_site", {}).get("peptide_chain", "A"),
    )
    print_stage(4, f"Inverse design ({method})")

    if method == "proteinmpnn":
        runner = SCRIPT_DIR / "run_mpnn_relax_all.py"
        method_arguments = [
            "--proteinmpnn-dir", str(paths["proteinmpnn_root"]),
            "--mpnn-env", config["envs"].get("proteinmpnn", "proteinmpnn"),
        ]
    elif method == "cyclicmpnn":
        runner = SCRIPT_DIR / "run_cyclicmpnn_all.py"
        method_arguments = [
            "--cyclicmpnn-dir", str(paths["cyclicmpnn_root"]),
            "--cyclicmpnn-env", config["envs"].get("cyclicmpnn", "cyclicmpnn"),
            "--model-name", design_config.get("model_name", "cyclicmpnn_48_010"),
        ]
    else:
        raise ValueError(
            f"Unknown inverse-design method {method!r}; "
            "expected 'proteinmpnn' or 'cyclicmpnn'."
        )

    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["biopython"],
        "python", str(runner),
        "--pdb-dir", str(paths["backbone_dir"]),
        "--output-dir", str(paths["mpnn_relax_dir"]),
        "--design-chain", design_chain,
        "--rounds", str(design_config["rounds"]),
        "--temperature", str(design_config.get("temperature", 0.1)),
        "--seed", str(design_config.get("seed", 1)),
        "--allowed-aas", design_config.get("allowed_aas", "CDEK"),
        "--pyrosetta-env", config["envs"].get("pyrosetta", "pyrosetta"),
        *method_arguments,
    ]
    if paths["distal_csv"] is not None:
        command.extend(["--scores", str(paths["distal_csv"])])
    if relax:
        command.extend(["--xml", str(paths["relax_xml"])])
    else:
        command.append("--no-relax")
    run_command(command)


def run_oracle(config, paths):
    print_stage(5, "AfCyc oracle")
    oracle_config = config["oracle"]
    command = [
        "conda", "run", "--no-capture-output",
        "-n", config["envs"]["biopython"],
        "python", str(SCRIPT_DIR / "run_afcyc_all.py"),
        "--design-dir", str(paths["mpnn_relax_dir"]),
        "--output-dir", str(paths["oracle_dir"]),
        "--target-pdb", str(paths["target_pdb"]),
        "--target-chain", oracle_config["target_chain"],
        "--params", str(paths["alphafold_params"]),
        "--afcyc-env", config["envs"]["afcyc"],
        "--recycles", str(oracle_config["recycles"]),
        "--seed", str(oracle_config["seed"]),
    ]
    env = os.environ.copy()
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    run_command(command, env=env)


def require_file(path, description):
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def require_directory(path, description):
    if not path.is_dir():
        raise FileNotFoundError(f"{description} not found: {path}")


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

    stage_settings = config.get("stages", {})
    stages = {
        "tensors": stage_settings.get("tensors", True),
        "rfdiffusion": stage_settings.get("rfdiffusion", True),
        "distal_site": stage_settings.get("distal_site", True),
        "inverse_design": stage_settings.get("inverse_design", True),
        "oracle": stage_settings.get("oracle", True),
    }

    design_config = config.get("mpnn_relax", {})
    method = design_config.get("method", "proteinmpnn")
    relax = design_config.get("relax", True)
    rounds = design_config.get("rounds")
    if stages["inverse_design"] and rounds is None:
        raise ValueError("[mpnn_relax].rounds is required")
    if stages["inverse_design"] and not relax and rounds != 1:
        raise ValueError("[mpnn_relax].relax = false requires rounds = 1")

    tensor_name = config.get("tensors", {}).get(
        "output_name",
        "interface_tensors_dir",
    )
    tensor_dir = project_dir / "tensors" / tensor_name
    generated_backbone_dir = project_dir / "rfdiffusion"
    generated_distal_csv = project_dir / "distal_site" / "distal_site.csv"

    input_config = config.get("inputs", {})
    if stages["rfdiffusion"]:
        backbone_dir = generated_backbone_dir
    else:
        backbone_dir = resolve_path(
            input_config.get("backbone_dir", generated_backbone_dir)
        )

    if stages["distal_site"]:
        distal_csv = generated_distal_csv
    elif "distal_scores" in input_config:
        distal_csv = resolve_path(input_config["distal_scores"])
    else:
        distal_csv = None

    paths = {
        "project_dir": project_dir,
        "target_pdb": target_pdb,
        "tensor_dir": tensor_dir,
        "tensor_list": tensor_dir.with_suffix(".txt"),
        "generated_backbone_dir": generated_backbone_dir,
        "backbone_dir": backbone_dir,
        "generated_distal_csv": generated_distal_csv,
        "distal_csv": distal_csv,
        "mpnn_relax_dir": project_dir / "mpnn_relax",
        "oracle_dir": project_dir / "oracle",
    }

    require_file(target_pdb, "Target PDB")
    if stages["rfdiffusion"]:
        paths["rfdiffusion_root"] = resolve_path(config["paths"]["rfdiffusion"])
        require_directory(paths["rfdiffusion_root"], "RFdiffusion directory")
    # Only validate backbones as an input when this run will not generate them.
    if (
        (stages["distal_site"] or stages["inverse_design"])
        and not stages["rfdiffusion"]
    ):
        require_directory(backbone_dir, "Backbone directory")

    # Likewise, an enabled distal-site stage will create its CSV before the
    # inverse-design stage needs it.
    if (
        stages["inverse_design"]
        and distal_csv is not None
        and not stages["distal_site"]
    ):
        require_file(distal_csv, "Distal-site CSV")

    if stages["inverse_design"] and method == "proteinmpnn":
        paths["proteinmpnn_root"] = resolve_path(config["paths"]["proteinmpnn"])
        require_directory(paths["proteinmpnn_root"], "ProteinMPNN directory")
    elif stages["inverse_design"] and method == "cyclicmpnn":
        paths["cyclicmpnn_root"] = resolve_path(config["paths"]["cyclicmpnn"])
        require_directory(paths["cyclicmpnn_root"], "CyclicMPNN directory")
    elif stages["inverse_design"]:
        raise ValueError(
            f"Unknown inverse-design method {method!r}; "
            "expected 'proteinmpnn' or 'cyclicmpnn'."
        )

    if stages["inverse_design"] and relax:
        paths["relax_xml"] = resolve_path(config["paths"]["relax_xml"])
        require_file(paths["relax_xml"], "Rosetta relax XML")
    if stages["oracle"]:
        paths["alphafold_params"] = resolve_path(config["paths"]["alphafold_params"])
        require_directory(paths["alphafold_params"], "AlphaFold parameters")

    output_directories = []
    if stages["tensors"]:
        output_directories.append(project_dir / "tensors")
    if stages["rfdiffusion"]:
        output_directories.append(generated_backbone_dir)
    if stages["distal_site"]:
        output_directories.append(generated_distal_csv.parent)
    if stages["inverse_design"]:
        output_directories.append(paths["mpnn_relax_dir"])
    if stages["oracle"]:
        output_directories.append(paths["oracle_dir"])
    for directory in output_directories:
        directory.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("Pipeline configuration")
    print("=" * 70)
    print(f"Config:        {config_path}")
    print(f"Project:       {project_dir}")
    print(f"Target PDB:    {target_pdb}")
    print(f"Backbones:     {backbone_dir}")
    print(f"Design method: {method}")
    print(f"Relax:         {relax}")
    for stage_name, enabled in stages.items():
        print(f"{stage_name:14} {'run' if enabled else 'skip'}")

    if stages["tensors"]:
        run_tensors(config, paths)
    if stages["rfdiffusion"]:
        run_rfdiffusion(config, paths)
    if stages["distal_site"]:
        run_distal_site(config, paths)
    if stages["inverse_design"]:
        run_inverse_design(config, paths)
    if stages["oracle"]:
        run_oracle(config, paths)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Project: {project_dir}")


if __name__ == "__main__":
    main()
