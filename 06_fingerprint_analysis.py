#!/usr/bin/env python3
"""
06_fingerprint_analysis.py

Analyze crystal fingerprints for generated crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes multiple crystal fingerprints
    * CrystalNNFingerprint
    * OPSiteFingerprint
    * VoronoiFingerprint
- Computes:
    * Mean fingerprint vector
    * Fingerprint norm
    * Pairwise cosine similarity
    * Nearest-neighbor similarity
    * PCA embedding
- Saves:
    * CSV with one row per structure
    * Individual fingerprint matrices
    * Similarity matrix
    * PCA coordinates
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF and PNG figures
- Prints dataset summary

Supported libraries
-------------------
- pymatgen
- matminer
- ase
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- seaborn
- tqdm
- joblib
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA

from pymatgen.core import Structure

from matminer.featurizers.site import (
    CrystalNNFingerprint,
    OPSiteFingerprint,
    VoronoiFingerprint,
)

from matminer.featurizers.structure import SiteStatsFingerprint

from common.results import (
    create_failed_result,
    create_structure_result,
)

from common.parallel import parallel_map

from common.statistics import (
    print_analysis_summary,
    print_top_categories,
    save_analysis_results,
)

from common.io import (
    create_analysis_output_directory,
    find_structure_files,
    load_structure,
    group_structure_files
)

from common.plotting import (
    configure_plotting,
    save_figure,
)

from common.cli import create_common_parser

configure_plotting()


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = create_common_parser(
        "Analyze crystal fingerprints."
    )

    return parser.parse_args()


def create_featurizers():
    """
    Create structure fingerprint featurizers.
    """

    crystalnn = SiteStatsFingerprint(
        CrystalNNFingerprint.from_preset("ops")
    )

    opsite = SiteStatsFingerprint(
        OPSiteFingerprint()
    )

    voronoi = SiteStatsFingerprint(
        VoronoiFingerprint()
    )

    return (
        crystalnn,
        opsite,
        voronoi,
    )


def featurize_structure(
    structure: Structure,
    featurizer,
):
    """
    Compute a structure fingerprint.

    Returns
    -------
    numpy.ndarray
    """

    return np.asarray(
        featurizer.featurize(
            structure,
        ),
        dtype=float,
    )


def compute_fingerprint_norm(
    fingerprint: np.ndarray,
):
    """
    Compute the Euclidean norm of a fingerprint.

    Parameters
    ----------
    fingerprint

    Returns
    -------
    float
    """

    return float(
        np.linalg.norm(
            fingerprint,
        )
    )


def compute_all_fingerprints(
    structure: Structure,
):
    """
    Compute all fingerprints for a structure.

    Returns
    -------
    dict
    """

    (
        crystalnn,
        opsite,
        voronoi,
    ) = create_featurizers()

    crystalnn_fp = featurize_structure(
        structure,
        crystalnn,
    )

    opsite_fp = featurize_structure(
        structure,
        opsite,
    )

    voronoi_fp = featurize_structure(
        structure,
        voronoi,
    )

    return {
        "crystalnn": crystalnn_fp,
        "opsite": opsite_fp,
        "voronoi": voronoi_fp,
    }


def analyze_structure(
    path: Path,
):
    """
    Analyze one crystal structure.

    Parameters
    ----------
    path

    Returns
    -------
    tuple
        analysis dictionary
        fingerprint dictionary
    """

    try:

        structure = load_structure(
            path,
        )

        fingerprints = compute_all_fingerprints(
            structure,
        )

        result = create_structure_result(
            path,
            structure,
            {
                "crystalnn_norm": compute_fingerprint_norm(
                    fingerprints["crystalnn"],
                ),
                "opsite_norm": compute_fingerprint_norm(
                    fingerprints["opsite"],
                ),
                "voronoi_norm": compute_fingerprint_norm(
                    fingerprints["voronoi"],
                ),
            },
        )

        return (
            result,
            fingerprints,
        )

    except Exception as exc:

        return (
            create_failed_result(
                path,
                {
                    "crystalnn_norm": np.nan,
                    "opsite_norm": np.nan,
                    "voronoi_norm": np.nan,
                },
                exc,
            ),
            None,
        )


def build_fingerprint_dataframe(
    fingerprints,
    key: str,
):
    """
    Convert one fingerprint type into a DataFrame.

    Parameters
    ----------
    fingerprints

    key

    Returns
    -------
    pandas.DataFrame
    """

    rows = []

    for file_name, fp_dict in fingerprints.items():

        vector = fp_dict[key]

        row = {
            "file": file_name,
        }

        row.update(
            {
                f"feature_{i:03d}": value
                for i, value in enumerate(vector)
            }
        )

        rows.append(
            row,
        )

    return pd.DataFrame(
        rows,
    )


def compute_similarity_statistics(
    matrix: pd.DataFrame,
):
    """
    Compute pairwise cosine similarities.

    Parameters
    ----------
    matrix

    Returns
    -------
    tuple
        similarity dataframe
        nearest-neighbor similarity
        mean similarity
    """

    feature_matrix = (
        matrix
        .drop(columns="file")
        .to_numpy(dtype=float)
    )

    similarity = cosine_similarity(
        feature_matrix,
    )

    similarity_df = pd.DataFrame(
        similarity,
        index=matrix["file"],
        columns=matrix["file"],
    )

    nearest_neighbor = []

    mean_similarity = []

    for row in similarity:

        values = row.copy()

        values = values[
            values < 0.999999
        ]

        if len(values):

            nearest_neighbor.append(
                float(
                    np.max(values)
                )
            )

            mean_similarity.append(
                float(
                    np.mean(values)
                )
            )

        else:

            nearest_neighbor.append(
                np.nan,
            )

            mean_similarity.append(
                np.nan,
            )

    return (
        similarity_df,
        nearest_neighbor,
        mean_similarity,
    )


def compute_pca(
    matrix: pd.DataFrame,
):
    """
    Compute PCA coordinates.

    Parameters
    ----------
    matrix

    Returns
    -------
    pandas.DataFrame
    """

    feature_matrix = (
        matrix
        .drop(columns="file")
        .to_numpy(dtype=float)
    )

    pca = PCA(
        n_components=2,
        random_state=42,
    )

    coordinates = pca.fit_transform(
        feature_matrix,
    )

    return pd.DataFrame(
        {
            "file": matrix["file"],
            "pc1": coordinates[:, 0],
            "pc2": coordinates[:, 1],
        }
    )


def print_summary(
    df: pd.DataFrame,
):
    """
    Print dataset summary.
    """

    if not print_analysis_summary(
        df,
        title="Fingerprint Analysis Summary",
        metrics=[
            "crystalnn_norm",
            "opsite_norm",
            "voronoi_norm",
            "mean_similarity",
            "nearest_neighbor_similarity",
        ],
        precision=4,
    ):
        return

    print_top_categories(
        df,
        column="reduced_formula",
        title="Most common formulas",
        top_n=10,
        width=20,
    )


def save_csv(
    analysis: pd.DataFrame,
    crystalnn: pd.DataFrame,
    opsite: pd.DataFrame,
    voronoi: pd.DataFrame,
    similarity: pd.DataFrame,
    pca: pd.DataFrame,
    output_dir: Path,
):
    """
    Save all CSV outputs.
    """

    save_analysis_results(
        analysis,
        output_dir,
        "fingerprint_analysis.csv",
        [
            "crystalnn_norm",
            "opsite_norm",
            "voronoi_norm",
            "mean_similarity",
            "nearest_neighbor_similarity",
        ],
    )

    crystalnn.to_csv(
        output_dir / "crystalnn_fingerprints.csv",
        index=False,
    )

    opsite.to_csv(
        output_dir / "opsite_fingerprints.csv",
        index=False,
    )

    voronoi.to_csv(
        output_dir / "voronoi_fingerprints.csv",
        index=False,
    )

    similarity.to_csv(
        output_dir / "similarity_matrix.csv",
    )

    pca.to_csv(
        output_dir / "pca_coordinates.csv",
        index=False,
    )


def plot_similarity_histogram(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot nearest-neighbor similarity histogram.
    """

    values = (
        df["nearest_neighbor_similarity"]
        .dropna()
    )

    if values.empty:
        return

    fig, ax = plt.subplots()

    sns.histplot(
        values,
        bins=40,
        kde=True,
        stat="count",
        ax=ax,
    )

    ax.set_xlabel(
        "Nearest Neighbor Cosine Similarity",
    )

    ax.set_ylabel(
        "Count",
    )

    ax.set_title(
        "Fingerprint Similarity Distribution",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_similarity_histogram",
    )


def plot_similarity_heatmap(
    similarity: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot fingerprint similarity heatmap.
    """

    if similarity.empty:
        return

    fig, ax = plt.subplots(
        figsize=(10, 8),
    )

    sns.heatmap(
        similarity,
        cmap="viridis",
        ax=ax,
        cbar_kws={
            "label": "Cosine Similarity",
        },
    )

    ax.set_xlabel(
        "Structure",
    )

    ax.set_ylabel(
        "Structure",
    )

    ax.set_title(
        "Fingerprint Similarity Matrix",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_heatmap",
    )


def plot_pca(
    pca: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot PCA embedding.
    """

    if pca.empty:
        return

    fig, ax = plt.subplots()

    ax.scatter(
        pca["pc1"],
        pca["pc2"],
        s=30,
    )

    ax.set_xlabel(
        "Principal Component 1",
    )

    ax.set_ylabel(
        "Principal Component 2",
    )

    ax.set_title(
        "Crystal Fingerprint PCA",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_pca",
    )


def plot_norm_boxplot(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot fingerprint norm distributions.
    """

    columns = [
        "crystalnn_norm",
        "opsite_norm",
        "voronoi_norm",
    ]

    plot_df = df[
        columns
    ].copy()

    if plot_df.dropna(
        how="all",
    ).empty:
        return

    fig, ax = plt.subplots()

    sns.boxplot(
        data=plot_df,
        ax=ax,
    )

    ax.set_xticklabels(
        [
            "CrystalNN",
            "OPSite",
            "Voronoi",
        ],
    )

    ax.set_ylabel(
        "Fingerprint Norm",
    )

    ax.set_title(
        "Fingerprint Norm Distribution",
    )

    save_figure(
        fig,
        output_dir,
        "fingerprint_boxplot",
    )


def create_all_plots(
    analysis: pd.DataFrame,
    similarity: pd.DataFrame,
    pca: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    valid = analysis[
        analysis["valid"]
    ]

    if valid.empty:
        return

    plot_similarity_histogram(
        valid,
        output_dir,
    )

    plot_similarity_heatmap(
        similarity,
        output_dir,
    )

    plot_pca(
        pca,
        output_dir,
    )

    plot_norm_boxplot(
        valid,
        output_dir,
    )



def main():

    args = parse_arguments()

    print()
    print("Searching for structures...")

    files = find_structure_files(
        args.input,
        args.recursive,
    )

    print(f"Found {len(files)} structures.")

    if not files:
        print("No supported files found.")
        return

    groups = group_structure_files(files)

    print(f"Found {len(groups)} compound groups.")

    for compound, compound_files in groups.items():

        print()
        print("=" * 60)
        print(f"Analyzing compound: {compound}")
        print(f"Structures: {len(compound_files)}")
        print("=" * 60)

        output_dir = create_analysis_output_directory(
            args.output / compound,
            "fingerprint",
        )

        results = parallel_map(
            compound_files,
            analyze_structure,
            args.workers,
            f"Computing {compound} fingerprints",
        )

        analysis_rows = []
        fingerprint_dict = {}

        for row, fingerprints in results:

            analysis_rows.append(row)

            if fingerprints is None:
                continue

            fingerprint_dict[row["file"]] = fingerprints

        analysis = pd.DataFrame(analysis_rows)

        valid = analysis[analysis["valid"]]

        if valid.empty:

            crystalnn = None
            opsite = None
            voronoi = None
            similarity = None
            pca = None

        else:

            crystalnn = build_fingerprint_dataframe(
                fingerprint_dict,
                "crystalnn",
            )

            opsite = build_fingerprint_dataframe(
                fingerprint_dict,
                "opsite",
            )

            voronoi = build_fingerprint_dataframe(
                fingerprint_dict,
                "voronoi",
            )

            (
                similarity,
                nearest_similarity,
                mean_similarity,
            ) = compute_similarity_statistics(
                crystalnn,
            )

            analysis.loc[
                analysis["valid"],
                "nearest_neighbor_similarity",
            ] = nearest_similarity

            analysis.loc[
                analysis["valid"],
                "mean_similarity",
            ] = mean_similarity

            pca = compute_pca(crystalnn)

            pca = pca.merge(
                analysis[
                    [
                        "file",
                        "formula",
                        "reduced_formula",
                    ]
                ],
                on="file",
                how="left",
            )

        print_summary(analysis)

        save_csv(
            analysis,
            crystalnn,
            opsite,
            voronoi,
            similarity,
            pca,
            output_dir,
        )

        print("Generating plots...")

        create_all_plots(
            analysis,
            similarity,
            pca,
            output_dir,
        )

        print()
        print(f"Results written to {output_dir.resolve()}")

    print()
    print("All compound groups analyzed.")



if __name__ == "__main__":
    main()
