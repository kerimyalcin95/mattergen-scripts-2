#!/usr/bin/env python3
"""
04_bond_analysis.py

Analyze bond lengths, bond angles and nearest-neighbor statistics.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes:
    * Bond length statistics
    * Bond angle statistics
    * Nearest-neighbor statistics
    * Mean coordination number
- Saves:
    * CSV with one row per structure
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF and PNG figures
- Prints summary to terminal

Supported libraries:
    pymatgen
    ase
    numpy
    pandas
    matplotlib
    seaborn
    scipy
    joblib
    tqdm
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from pymatgen.core import Structure
from pymatgen.core.local_env import CrystalNN

from common.results import (
    create_failed_result,
    create_structure_result,
)

from common.parallel import parallel_map

from common.io import (
    create_analysis_output_directory,
    find_structure_files,
    load_structure,
)

from common.statistics import (
    print_analysis_summary,
    print_top_categories,
    save_analysis_results,
    statistics_to_dict,
)

from common.cli import create_common_parser

from common.plotting import (
    plot_histogram,
    save_figure,
    configure_plotting
)

configure_plotting()


def main():

    args = parse_arguments()

    output_dir = create_analysis_output_directory(
        args.output,
        "bond",
    )

    print()
    print("Searching for structures...")

    files = find_structure_files(
        args.input,
        args.recursive,
    )

    print(f"Found {len(files)} structures.")

    if len(files) == 0:
        print("No supported files found.")
        return

    print()

    rows = parallel_map(
        files,
        analyze_structure,
        args.workers,
        "Analyzing bonds",
    )

    df = pd.DataFrame(
        rows,
    )

    print_summary(
        df,
    )

    save_csv(
        df,
        output_dir,
    )

    print("Creating figures...")

    create_all_plots(
        df,
        output_dir,
    )

    print()
    print(f"Results written to {output_dir.resolve()}")
    print("Done.")


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = create_common_parser(
        "Analyze bond lengths, bond angles and nearest-neighbor statistics."
    )

    return parser.parse_args()


def unique_pair(
    index_a: int,
    index_b: int,
):
    """
    Return an order-independent atom pair.
    """

    return tuple(
        sorted(
            (
                index_a,
                index_b,
            )
        )
    )


def extract_bonds(
    structure: Structure,
):
    """
    Extract unique bonds using CrystalNN.

    Returns
    -------
    tuple[list[float], list[int], dict[int, list[int]]]
        bond_lengths
        coordination_numbers
        neighbor_map
    """

    crystal_nn = CrystalNN(
        distance_cutoffs=None,
        x_diff_weight=0,
        porous_adjustment=False,
    )

    unique_bonds = set()

    bond_lengths = []

    coordination_numbers = []

    neighbor_map = {
        index: []
        for index in range(len(structure))
    }

    for site_index in range(len(structure)):

        try:

            neighbors = crystal_nn.get_nn_info(
                structure,
                site_index,
            )

        except Exception:

            coordination_numbers.append(np.nan)

            continue

        coordination_numbers.append(
            len(neighbors)
        )

        for neighbor in neighbors:

            neighbor_index = neighbor["site_index"]

            pair = unique_pair(
                site_index,
                neighbor_index,
            )

            neighbor_map[site_index].append(
                neighbor_index
            )

            if pair in unique_bonds:
                continue

            unique_bonds.add(pair)

            distance = structure.get_distance(
                site_index,
                neighbor_index,
            )

            bond_lengths.append(
                float(distance)
            )

    return (
        bond_lengths,
        coordination_numbers,
        neighbor_map,
    )


def compute_nearest_neighbor_distances(
    structure: Structure,
    neighbor_map,
):
    """
    Compute nearest-neighbor distance for every atom.

    Returns
    -------
    list[float]
    """

    nearest_distances = []

    for atom_index in range(len(structure)):

        neighbors = neighbor_map.get(
            atom_index,
            [],
        )

        if len(neighbors) == 0:

            nearest_distances.append(
                np.nan,
            )

            continue

        distances = [
            structure.get_distance(
                atom_index,
                neighbor_index,
            )
            for neighbor_index in neighbors
        ]

        nearest_distances.append(
            float(min(distances))
        )

    return nearest_distances


def summarize_bond_statistics(
    bond_lengths,
    nearest_neighbor_distances,
    coordination_numbers,
):
    """
    Compute summary statistics.

    Returns
    -------
    dict
    """

    results = {}

    results.update(
        statistics_to_dict(
            bond_lengths,
            "bond",
        )
    )

    results.update(
        statistics_to_dict(
            nearest_neighbor_distances,
            "nearest_neighbor",
        )
    )

    results.update(
        statistics_to_dict(
            coordination_numbers,
            "coordination",
        )
    )

    results["number_of_bonds"] = len(
        bond_lengths,
    )

    return results


def compute_bond_angles(
    structure: Structure,
    neighbor_map,
):
    """
    Compute all bond angles in the structure.

    For every central atom B and neighboring atoms
    A and C, the angle A-B-C is calculated.

    Returns
    -------
    list[float]
        Bond angles in degrees.
    """

    bond_angles = []

    for center_index in range(len(structure)):

        neighbors = neighbor_map.get(
            center_index,
            [],
        )

        if len(neighbors) < 2:
            continue

        center = structure[center_index].coords

        for i in range(len(neighbors) - 1):

            first_index = neighbors[i]

            vector_1 = (
                structure[first_index].coords
                - center
            )

            norm_1 = np.linalg.norm(
                vector_1,
            )

            if norm_1 == 0:
                continue

            for j in range(i + 1, len(neighbors)):

                second_index = neighbors[j]

                vector_2 = (
                    structure[second_index].coords
                    - center
                )

                norm_2 = np.linalg.norm(
                    vector_2,
                )

                if norm_2 == 0:
                    continue

                cosine = np.dot(
                    vector_1,
                    vector_2,
                ) / (
                    norm_1 * norm_2
                )

                cosine = np.clip(
                    cosine,
                    -1.0,
                    1.0,
                )

                angle = np.degrees(
                    np.arccos(
                        cosine,
                    )
                )

                bond_angles.append(
                    float(angle)
                )

    return bond_angles


def summarize_angle_statistics(
    bond_angles,
):
    """
    Compute bond-angle statistics.

    Returns
    -------
    dict
    """

    results = statistics_to_dict(
        bond_angles,
        "bond_angle",
    )

    results["number_of_angles"] = len(
        bond_angles,
    )

    return results


def collect_bond_information(
    structure: Structure,
):
    """
    Compute all bond-related information for one structure.

    Returns
    -------
    dict
    """

    (
        bond_lengths,
        coordination_numbers,
        neighbor_map,
    ) = extract_bonds(
        structure,
    )

    nearest_neighbor_distances = (
        compute_nearest_neighbor_distances(
            structure,
            neighbor_map,
        )
    )

    bond_angles = compute_bond_angles(
        structure,
        neighbor_map,
    )

    results = summarize_bond_statistics(
        bond_lengths,
        nearest_neighbor_distances,
        coordination_numbers,
    )

    results.update(
        summarize_angle_statistics(
            bond_angles,
        )
    )

    return results


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
    dict
    """

    try:

        structure = load_structure(
            path,
        )

        statistics = collect_bond_information(
            structure,
        )

        result = create_structure_result(
            path,
            structure,
            statistics,
        )

        return result

    except Exception as exc:

        return create_failed_result(
            path,
            {
                "bond_mean": np.nan,
                "bond_median": np.nan,
                "bond_std": np.nan,
                "bond_min": np.nan,
                "bond_max": np.nan,
                "nearest_neighbor_mean": np.nan,
                "nearest_neighbor_median": np.nan,
                "nearest_neighbor_std": np.nan,
                "nearest_neighbor_min": np.nan,
                "nearest_neighbor_max": np.nan,
                "coordination_mean": np.nan,
                "coordination_median": np.nan,
                "coordination_std": np.nan,
                "coordination_min": np.nan,
                "coordination_max": np.nan,
                "bond_angle_mean": np.nan,
                "bond_angle_median": np.nan,
                "bond_angle_std": np.nan,
                "bond_angle_min": np.nan,
                "bond_angle_max": np.nan,
                "number_of_bonds": np.nan,
                "number_of_angles": np.nan,
            },
            exc,
        )


def print_summary(
    df: pd.DataFrame,
):
    """
    Print a dataset summary.
    """

    if not print_analysis_summary(
        df,
        title="Bond Analysis Summary",
        metrics=[
            "bond_mean",
            "bond_angle_mean",
            "nearest_neighbor_mean",
            "coordination_mean",
        ],
        precision=3,
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
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Save analysis results.
    """

    save_analysis_results(
        df,
        output_dir,
        "bond_analysis.csv",
        [
            "bond_mean",
            "bond_angle_mean",
            "nearest_neighbor_mean",
            "coordination_mean",
            "number_of_bonds",
            "number_of_angles",
        ],
    )


def plot_boxplots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot boxplots of the main statistics.
    """

    columns = [
        "bond_mean",
        "bond_angle_mean",
        "nearest_neighbor_mean",
        "coordination_mean",
    ]

    labels = [
        "Bond\nLength",
        "Bond\nAngle",
        "Nearest\nNeighbor",
        "Coordination",
    ]

    plot_df = df[
        columns
    ].copy()

    if plot_df.dropna(how="all").empty:
        return

    fig, ax = plt.subplots(
        figsize=(8, 5),
    )

    sns.boxplot(
        data=plot_df,
        ax=ax,
    )

    ax.set_xticklabels(
        labels,
    )

    ax.set_ylabel("Value")

    ax.set_title(
        "Bond Statistics",
    )

    save_figure(
        fig,
        output_dir,
        "bond_statistics_boxplot",
    )


def create_all_plots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Create every figure.
    """

    valid = df[
        df["valid"]
    ]

    if valid.empty:
        return

    plot_histogram(
        valid["bond_mean"],
        xlabel="Mean Bond Length (Å)",
        title="Bond Length Distribution",
        filename="bond_length_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["bond_angle_mean"],
        xlabel="Mean Bond Angle (°)",
        title="Bond Angle Distribution",
        filename="bond_angle_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["nearest_neighbor_mean"],
        xlabel="Nearest Neighbor Distance (Å)",
        title="Nearest Neighbor Distribution",
        filename="nearest_neighbor_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["coordination_mean"],
        xlabel="Coordination Number",
        title="Coordination Number Distribution",
        filename="coordination_histogram",
        output_dir=output_dir,
        bins=20,
    )

    plot_boxplots(
        valid,
        output_dir,
    )


if __name__ == "__main__":
    main()
