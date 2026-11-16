#!/usr/bin/env python3
"""
01_basic_statistics.py

Compute basic statistics for a dataset of crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes:
    * Chemical formula
    * Number of atoms
    * Lattice parameters
    * Cell volume
    * Density
    * Space group
- Saves:
    * CSV with one row per structure
    * Summary CSV
    * Publication-quality PDF figures
- Prints summary to terminal

Supported libraries:
    pymatgen
    ase
    numpy
    pandas
    scipy
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

from pymatgen.core import Structure

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from common.parallel import parallel_map

from common.results import (
    create_failed_result,
    create_structure_result,
)

from common.io import (
    ensure_output_directory,
    find_structure_files,
    load_structure,
)

from common.statistics import (
    print_analysis_summary,
    print_top_categories,
    save_analysis_results,
)

from common.cli import create_common_parser

from common.plotting import (
    configure_plotting,
    plot_category_distribution,
    plot_histogram,
)

configure_plotting()


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = create_common_parser(
        "Compute basic statistics of crystal datasets."
    )

    return parser.parse_args()


def safe_spacegroup(
    structure: Structure,
):
    """
    Determine the space group.

    Returns
    -------
    str
    """

    try:
        analyzer = SpacegroupAnalyzer(
            structure,
            symprec=0.1,
        )

        return analyzer.get_space_group_symbol()

    except Exception:
        return "Unknown"


def safe_density(
    structure: Structure,
):
    """
    Compute density.

    Returns
    -------
    float
    """

    try:
        return float(structure.density)

    except Exception:
        return np.nan


def safe_volume(
    structure: Structure,
):
    """
    Compute cell volume.
    """

    try:
        return float(structure.volume)

    except Exception:
        return np.nan


def analyze_structure(path: Path):
    """
    Analyze a single structure.

    Parameters
    ----------
    path : Path

    Returns
    -------
    dict
        Dictionary containing all computed properties.
    """

    try:
        structure = load_structure(path)

        lattice = structure.lattice

        result = create_structure_result(
            path,
            structure,
            {
                "num_elements": len(structure.composition.elements),
                "a": lattice.a,
                "b": lattice.b,
                "c": lattice.c,
                "alpha": lattice.alpha,
                "beta": lattice.beta,
                "gamma": lattice.gamma,
                "volume": safe_volume(structure),
                "density": safe_density(structure),
                "space_group": safe_spacegroup(structure),
            },
        )

        return result

    except Exception as exc:

        return create_failed_result(
            path,
            {
                "num_elements": np.nan,
                "a": np.nan,
                "b": np.nan,
                "c": np.nan,
                "alpha": np.nan,
                "beta": np.nan,
                "gamma": np.nan,
                "volume": np.nan,
                "density": np.nan,
                "space_group": "Invalid",
            },
            exc,
        )


def print_summary(
    df: pd.DataFrame,
):
    """
    Print a dataset summary.
    """

    print_analysis_summary(
        df,
        title="MatterGen Dataset Summary",
        metrics=[
            "num_atoms",
            "volume",
            "density",
            "a",
            "b",
            "c",
        ],
        precision=3,
    )

    print_top_categories(
        df,
        column="space_group",
        title="Most common space groups",
        top_n=10,
        width=12,
    )


def save_csv(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Save the per-structure statistics.
    """

    save_analysis_results(
        df,
        output_dir,
        "basic_statistics.csv",
        [
            "num_atoms",
            "volume",
            "density",
            "a",
            "b",
            "c",
            "alpha",
            "beta",
            "gamma",
        ],
    )


def plot_lattice_parameters(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot a, b and c distributions.
    """

    plt.figure(figsize=(10, 6))

    for column in ["a", "b", "c"]:

        sns.histplot(
            df[column].dropna(),
            bins=30,
            kde=True,
            label=column,
            stat="density",
            element="step",
            fill=False,
        )

    plt.xlabel("Lattice parameter (Å)")
    plt.ylabel("Density")
    plt.title("Lattice parameter distributions")

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_dir / "lattice_parameters.pdf",
        bbox_inches="tight",
    )

    plt.savefig(
        output_dir / "lattice_parameters.png",
        bbox_inches="tight",
    )

    plt.close()


def create_all_plots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    plot_histogram(
        df["num_atoms"],
        xlabel="Atoms per cell",
        title="Atoms per cell",
        filename="num_atoms_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        df["volume"],
        xlabel="Cell volume (Å³)",
        title="Cell volume (Å³)",
        filename="volume_histogram",
        output_dir=output_dir,
    )

    plot_histogram(
        df["density"],
        xlabel="Density (g/cm³)",
        title="Density (g/cm³)",
        filename="density_histogram",
        output_dir=output_dir,
    )

    plot_lattice_parameters(
        df,
        output_dir,
    )

    plot_category_distribution(
        df["space_group"],
        xlabel="Count",
        ylabel="Space group",
        title="Most common space groups",
        filename="spacegroup_distribution",
        output_dir=output_dir,
        top_n=20,
    )

    plot_category_distribution(
        df["reduced_formula"],
        xlabel="Count",
        ylabel="Formula",
        title="Most common compositions",
        filename="composition_distribution",
        output_dir=output_dir,
        top_n=20,
    )


def main():

    args = parse_arguments()

    ensure_output_directory(
        args.output,
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
        "Analyzing structures",
    )

    df = pd.DataFrame(
        rows,
    )

    print_summary(df)

    save_csv(
        df,
        args.output,
    )

    print("Generating plots...")

    create_all_plots(
        df,
        args.output,
    )

    print()

    print(f"Results written to {args.output.resolve()}")

    print("Done.")


if __name__ == "__main__":
    main()
