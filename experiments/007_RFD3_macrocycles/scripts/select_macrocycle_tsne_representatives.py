#!/usr/bin/env python3

"""Select one high-quality AfCyc prediction from each t-SNE/GMM cluster.

A prediction qualifies only when its backbone RMSD is below 1 Angstrom and
its AfCyc pLDDT is above 0.8.  Within each cluster, the qualifying prediction
with the highest pLDDT is selected; lower RMSD breaks a tie.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PLDDT_CUTOFF = 0.8
RMSD_CUTOFF = 1.0
EXPECTED_CLUSTERS = 40

SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1] / "monomer_self_consistency"
DEFAULT_STRUCTURAL_TSNE_DIR = EXPERIMENT_DIR / "structural_tsne"
DEFAULT_ORACLE_ANALYSIS_DIR = EXPERIMENT_DIR / "analysis" / "run_01"
DEFAULT_OUTPUT_DIR = DEFAULT_STRUCTURAL_TSNE_DIR / "pymol_examples" / "run_01"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structural-tsne-dir",
        type=Path,
        default=DEFAULT_STRUCTURAL_TSNE_DIR,
        help=f"Structural t-SNE directory. Default: {DEFAULT_STRUCTURAL_TSNE_DIR}",
    )
    parser.add_argument(
        "--oracle-analysis-dir",
        type=Path,
        default=DEFAULT_ORACLE_ANALYSIS_DIR,
        help=f"AfCyc analysis directory. Default: {DEFAULT_ORACLE_ANALYSIS_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Manifest output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    structural_tsne_dir = args.structural_tsne_dir.expanduser().resolve()
    oracle_analysis_dir = args.oracle_analysis_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    clusters = pd.read_csv(
        structural_tsne_dir / "embedding" / "tsne_clusters.tsv",
        sep="\t",
    )
    predictions = pd.read_csv(
        oracle_analysis_dir / "prediction_results.tsv",
        sep="\t",
    )

    qualifying = predictions[
        (predictions["plddt"] > PLDDT_CUTOFF)
        & (predictions["backbone_rmsd_angstrom"] < RMSD_CUTOFF)
    ].copy()

    qualifying = qualifying.merge(
        clusters[
            [
                "backbone",
                "peptide_length",
                "cluster",
                "tsne_component_1",
                "tsne_component_2",
            ]
        ],
        on=["backbone", "peptide_length"],
        how="inner",
        validate="many_to_one",
    )

    qualifying = qualifying.sort_values(
        by=[
            "cluster",
            "plddt",
            "backbone_rmsd_angstrom",
            "backbone",
            "sequence_id",
        ],
        ascending=[True, False, True, True, True],
    )

    selected_by_cluster = {
        int(cluster): rows.iloc[0]
        for cluster, rows in qualifying.groupby("cluster", sort=True)
    }

    manifest_rows = []
    for cluster in range(1, EXPECTED_CLUSTERS + 1):
        if cluster not in selected_by_cluster:
            manifest_rows.append(
                {
                    "cluster": cluster,
                    "status": "no_qualifying_pair",
                    "backbone": "",
                    "peptide_length": "",
                    "sequence_id": "",
                    "sequence": "",
                    "pLDDT": "",
                    "backbone_RMSD": "",
                    "design_path": "",
                    "AfCycDesign_prediction_path": "",
                    "tsne_component_1": "",
                    "tsne_component_2": "",
                }
            )
            continue

        row = selected_by_cluster[cluster]
        manifest_rows.append(
            {
                "cluster": cluster,
                "status": "selected",
                "backbone": row["backbone"],
                "peptide_length": int(row["peptide_length"]),
                "sequence_id": int(row["sequence_id"]),
                "sequence": row["sequence"],
                "pLDDT": row["plddt"],
                "backbone_RMSD": row["backbone_rmsd_angstrom"],
                "design_path": row["source_path"],
                "AfCycDesign_prediction_path": row["prediction_path"],
                "tsne_component_1": row["tsne_component_1"],
                "tsne_component_2": row["tsne_component_2"],
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "cluster_representatives.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    missing_clusters = manifest.loc[
        manifest["status"] == "no_qualifying_pair",
        "cluster",
    ].tolist()

    print(f"Qualifying predictions: {len(qualifying):,}")
    print(f"Clusters with a representative: {EXPECTED_CLUSTERS - len(missing_clusters)}")
    if missing_clusters:
        print(f"Clusters without a qualifying pair: {missing_clusters}")
    else:
        print("Every cluster has a qualifying representative.")
    print(f"Saved: {manifest_path}")


if __name__ == "__main__":
    main()
