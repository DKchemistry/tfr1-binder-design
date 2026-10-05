#!/usr/bin/env python3

"""Embed the macrocycle TM-score matrix with t-SNE and cluster it with GMM.

Each row of the 1,200 by 1,200 TM-score matrix is treated directly as the
high-dimensional representation of one backbone.  The matrix is not converted
to a precomputed distance matrix.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.manifold import TSNE
from sklearn.mixture import GaussianMixture


SCRIPT_PATH = Path(__file__).resolve()
EXPERIMENT_DIR = SCRIPT_PATH.parents[1] / "monomer_self_consistency"
DEFAULT_ANALYSIS_DIR = EXPERIMENT_DIR / "structural_tsne"
DEFAULT_RANDOM_SEED = 0
GMM_COMPONENTS = 40


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=DEFAULT_ANALYSIS_DIR,
        help=f"Structural t-SNE directory. Default: {DEFAULT_ANALYSIS_DIR}",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=DEFAULT_RANDOM_SEED,
        help=f"Random seed for t-SNE and GMM. Default: {DEFAULT_RANDOM_SEED}",
    )
    return parser.parse_args()


def set_plot_style() -> None:
    """Use the compact academic style of the other experiment 007 plots."""

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def make_cluster_colors() -> list[tuple[float, float, float, float]]:
    """Return 40 categorical colors from Matplotlib palettes."""

    colors = list(plt.get_cmap("tab20").colors)
    colors.extend(plt.get_cmap("tab20b").colors)
    return colors


def plot_embedding(results: pd.DataFrame, output_dir: Path) -> None:
    """Plot t-SNE coordinates colored by the 40 GMM clusters."""

    set_plot_style()
    colors = make_cluster_colors()
    fig, ax = plt.subplots(figsize=(5.2, 4.4))

    for cluster in range(1, GMM_COMPONENTS + 1):
        cluster_rows = results[results["cluster"] == cluster]
        ax.scatter(
            cluster_rows["tsne_component_1"],
            cluster_rows["tsne_component_2"],
            s=13,
            alpha=0.8,
            edgecolors="none",
            color=colors[cluster - 1],
        )

        # Cluster numbers make the plot usable with the inspection manifest
        # without requiring a large 40-entry legend.
        ax.text(
            cluster_rows["tsne_component_1"].median(),
            cluster_rows["tsne_component_2"].median(),
            str(cluster),
            ha="center",
            va="center",
            fontsize=6,
            weight="bold",
            bbox={
                "boxstyle": "circle,pad=0.15",
                "facecolor": "white",
                "edgecolor": colors[cluster - 1],
                "linewidth": 0.7,
                "alpha": 0.9,
            },
        )

    ax.set_xlabel("t-SNE component 1")
    ax.set_ylabel("t-SNE component 2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3)

    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "macrocycle_structural_tsne.pdf"
    png_path = output_dir / "macrocycle_structural_tsne.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


def json_safe_parameters(parameters: dict[str, object]) -> dict[str, object]:
    """Convert estimator parameters to values accepted by json.dump()."""

    safe_parameters = {}
    for key, value in parameters.items():
        if isinstance(value, np.generic):
            value = value.item()
        safe_parameters[key] = value
    return safe_parameters


def main() -> None:
    args = parse_arguments()
    analysis_dir = args.analysis_dir.expanduser().resolve()

    matrix_path = analysis_dir / "tmalign" / "tm_scores.npy"
    index_path = analysis_dir / "structure_index.tsv"

    matrix = np.load(matrix_path)
    structure_index = pd.read_csv(index_path, sep="\t")

    if matrix.shape != (len(structure_index), len(structure_index)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match "
            f"{len(structure_index)} indexed structures"
        )

    tsne = TSNE(
        n_components=2,
        random_state=args.random_seed,
    )
    coordinates = tsne.fit_transform(matrix)

    gmm = GaussianMixture(
        n_components=GMM_COMPONENTS,
        random_state=args.random_seed,
    )
    zero_based_clusters = gmm.fit_predict(coordinates)
    membership_probabilities = gmm.predict_proba(coordinates).max(axis=1)

    results = structure_index.copy()
    results["tsne_component_1"] = coordinates[:, 0]
    results["tsne_component_2"] = coordinates[:, 1]
    results["cluster"] = zero_based_clusters + 1
    results["cluster_probability"] = membership_probabilities

    embedding_dir = analysis_dir / "embedding"
    embedding_dir.mkdir(parents=True, exist_ok=True)
    results_path = embedding_dir / "tsne_clusters.tsv"
    results.to_csv(results_path, sep="\t", index=False)

    metadata = {
        "input": "rows of the symmetric averaged TM-score matrix",
        "scikit_learn_version": sklearn.__version__,
        "tsne_parameters": json_safe_parameters(tsne.get_params()),
        "gmm_parameters": json_safe_parameters(gmm.get_params()),
    }
    metadata_path = embedding_dir / "metadata.json"
    with metadata_path.open("w") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")

    print("Cluster sizes:")
    cluster_sizes = results["cluster"].value_counts().sort_index()
    for cluster, size in cluster_sizes.items():
        print(f"  Cluster {cluster:>2}: {size:>3} backbones")

    print(f"Saved: {results_path}")
    print(f"Saved: {metadata_path}")
    plot_embedding(results, analysis_dir / "figures")


if __name__ == "__main__":
    main()
