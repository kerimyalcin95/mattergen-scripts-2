#!/usr/bin/env python3
"""
03_coordination_analysis.py

Analyze local coordination environments of crystal structures using both
CrystalNN and VoronoiNN.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes per-atom coordination numbers with:
    * CrystalNN
    * VoronoiNN
- Computes per-structure statistics:
    * mean, median, std, min, max coordination
    * agreement fraction between both methods
    * mean and max absolute difference
- Saves:
    * CSV with one row per structure
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF/PNG figures
- Prints a human-readable summary to terminal

Supported libraries:
    pymatgen
    ase
    numpy
    pandas
    matplotlib
    seaborn
    joblib
    tqdm
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from pymatgen.core.local_env import CrystalNN, VoronoiNN

from common.parallel import parallel_map

from common.results import (
    create_failed_result,
    create_structure_result,
)

from common.io import (
    create_analysis_output_directory,
    find_structure_files,
    load_structure,
    group_structure_files
)

from common.statistics import (
    print_analysis_summary,
    print_top_categories,
    save_analysis_results,
    statistics_to_dict,
)

from common.cli import create_common_parser

from common.plotting import (
    plot_category_distribution,
    plot_histogram,
    save_figure,
    configure_plotting
)

configure_plotting()


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = create_common_parser(
        "Analyze coordination environments using CrystalNN and VoronoiNN."
    )

    return parser.parse_args()


def coordination_number_from_nninfo(nn_info):
    """
    Convert neighbor information into an integer coordination number.
    """

    if nn_info is None:
        return 0
    return len(nn_info)


def analyze_structure(path: Path):
    """
    Analyze coordination for one structure.
    """

    try:
        structure = load_structure(path)
        crystal_nn = CrystalNN(
            distance_cutoffs=None,
            x_diff_weight=0,
            porous_adjustment=False,
        )
        voronoi_nn = VoronoiNN()

        crystal_cns = []
        voronoi_cns = []
        per_site_differences = []
        agreement_hits = 0

        for site_index in range(len(structure)):
            try:
                c_nn = crystal_nn.get_nn_info(structure, site_index)
                c_cn = coordination_number_from_nninfo(c_nn)
            except Exception:
                c_cn = np.nan

            try:
                v_nn = voronoi_nn.get_nn_info(structure, site_index)
                v_cn = coordination_number_from_nninfo(v_nn)
            except Exception:
                v_cn = np.nan

            crystal_cns.append(c_cn)
            voronoi_cns.append(v_cn)

            if np.isfinite(c_cn) and np.isfinite(v_cn):
                diff = abs(float(c_cn) - float(v_cn))
                per_site_differences.append(diff)
                if diff == 0:
                    agreement_hits += 1

        crystal_values = pd.Series(crystal_cns, dtype="float64").dropna()
        voronoi_values = pd.Series(voronoi_cns, dtype="float64").dropna()
        diff_values = pd.Series(per_site_differences, dtype="float64").dropna()

        formula = structure.composition.formula
        reduced_formula = structure.composition.reduced_formula

        crystalnn_statistics = statistics_to_dict(
            crystal_values,
            "crystalnn",
            "_coordination",
        )

        voronoi_statistics = statistics_to_dict(
            voronoi_values,
            "voronoi",
            "_coordination",
        )

        statistics = {}

        statistics.update(
            crystalnn_statistics,
        )

        statistics.update(
            voronoi_statistics,
        )

        statistics.update(
            {
                "mean_absolute_difference": (
                    diff_values.mean()
                    if not diff_values.empty
                    else np.nan
                ),
                "max_absolute_difference": (
                    diff_values.max()
                    if not diff_values.empty
                    else np.nan
                ),
                "agreement_fraction": (
                    agreement_hits / len(structure)
                    if len(structure) > 0
                    else np.nan
                ),
            }
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
                "crystalnn_mean_coordination": np.nan,
                "crystalnn_median_coordination": np.nan,
                "crystalnn_std_coordination": np.nan,
                "crystalnn_min_coordination": np.nan,
                "crystalnn_max_coordination": np.nan,
                "voronoi_mean_coordination": np.nan,
                "voronoi_median_coordination": np.nan,
                "voronoi_std_coordination": np.nan,
                "voronoi_min_coordination": np.nan,
                "voronoi_max_coordination": np.nan,
                "mean_absolute_difference": np.nan,
                "max_absolute_difference": np.nan,
                "agreement_fraction": np.nan,
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
        title="Coordination Analysis Summary",
        metrics=[
            "crystalnn_mean_coordination",
            "voronoi_mean_coordination",
            "mean_absolute_difference",
            "agreement_fraction",
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
        "coordination_analysis.csv",
        [
            "crystalnn_mean_coordination",
            "voronoi_mean_coordination",
            "mean_absolute_difference",
            "agreement_fraction",
        ],
    )


def plot_comparison(df: pd.DataFrame, output_dir: Path):
    """
    Scatter plot comparing CrystalNN and VoronoiNN.
    """

    valid = df[df["valid"]].dropna(
        subset=["crystalnn_mean_coordination", "voronoi_mean_coordination"])
    if valid.empty:
        return

    plt.figure(figsize=(6, 6))
    sns.scatterplot(
        data=valid,
        x="crystalnn_mean_coordination",
        y="voronoi_mean_coordination",
        s=40,
    )
    max_value = max(
        valid["crystalnn_mean_coordination"].max(),
        valid["voronoi_mean_coordination"].max(),
    )
    plt.plot([0, max_value], [0, max_value], linestyle="--")
    plt.xlabel("CrystalNN mean coordination")
    plt.ylabel("VoronoiNN mean coordination")
    plt.title("CrystalNN vs VoronoiNN")

    save_figure(
        plt.gcf(),
        output_dir,
        "coordination_scatter",
    )


def plot_difference(df: pd.DataFrame, output_dir: Path):
    """
    Plot distribution of coordination differences.
    """

    values = df.loc[df["valid"], "mean_absolute_difference"].dropna()
    if values.empty:
        return

    plt.figure()
    sns.histplot(values, bins=30, kde=True)
    plt.xlabel("Mean absolute difference")
    plt.ylabel("Count")
    plt.title("CrystalNN vs VoronoiNN difference")

    plt.tight_layout()

    for ext in ("pdf", "png"):
        plt.savefig(
            output_dir /
            f"coordination_difference_histogram.{ext}",
            bbox_inches="tight")
    plt.close()


def plot_boxplot(df: pd.DataFrame, output_dir: Path):
    """
    Compare coordination number distributions.
    """

    valid = df[df["valid"]]
    if valid.empty:
        return

    melted = pd.melt(
        valid,
        value_vars=[
            "crystalnn_mean_coordination",
            "voronoi_mean_coordination"],
        var_name="method",
        value_name="coordination",
    )

    plt.figure(figsize=(7, 5))
    sns.boxplot(data=melted, x="method", y="coordination")
    plt.xlabel("Method")
    plt.ylabel("Mean coordination")
    plt.title("Coordination comparison")

    save_figure(
        plt.gcf(),
        output_dir,
        "coordination_boxplot",
    )


def create_all_plots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    plot_histogram(
        df["crystalnn_mean_coordination"],
        xlabel="CrystalNN mean coordination",
        title="CrystalNN mean coordination",
        filename="crystalnn_coordination_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        df["voronoi_mean_coordination"],
        xlabel="VoronoiNN mean coordination",
        title="VoronoiNN mean coordination",
        filename="voronoi_coordination_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        df["agreement_fraction"],
        xlabel="Agreement fraction",
        title="Agreement fraction",
        filename="agreement_fraction_histogram",
        output_dir=output_dir,
    )

    plot_difference(
        df,
        output_dir,
    )

    plot_comparison(
        df,
        output_dir,
    )

    plot_boxplot(
        df,
        output_dir,
    )

    plot_category_distribution(
        df.loc[
            df["valid"],
            "reduced_formula",
        ],
        xlabel="Count",
        ylabel="Formula",
        title="Most common compositions",
        filename="composition_distribution",
        output_dir=output_dir,
        top_n=20,
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
            "coordination",
        )

        rows = parallel_map(
            compound_files,
            analyze_structure,
            args.workers,
            f"Analyzing {compound} coordination",
        )

        df = pd.DataFrame(rows)

        print_summary(df)

        save_csv(
            df,
            output_dir,
        )

        print("Generating plots...")

        create_all_plots(
            df,
            output_dir,
        )

        print()
        print(f"Results written to {output_dir.resolve()}")

    print()
    print("All compound groups analyzed.")



if __name__ == "__main__":
    main()
