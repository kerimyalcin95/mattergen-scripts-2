#!/usr/bin/env python3
"""
02_symmetry_analysis.py

Analyze the crystallographic symmetry of a dataset of crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes:
    * Space group number
    * Space group symbol
    * Crystal system
    * Bravais lattice
    * Point group
    * Hall symbol
    * Number of symmetry operations
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
    matplotlib
    seaborn
    joblib
    tqdm
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from functools import partial

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from common.parallel import parallel_map

from common.io import (
    create_analysis_output_directory,
    find_structure_files,
    load_structure,
    group_structure_files,

)

from common.statistics import (
    print_analysis_summary,
    print_top_categories,
    save_analysis_results,
)

from common.cli import create_common_parser

from common.plotting import (
    plot_category_distribution,
    plot_histogram,
    configure_plotting
)

configure_plotting()

DEFAULT_SYMPREC = 0.1


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = create_common_parser(
        "Analyze crystallographic symmetry."
    )

    parser.add_argument(
        "--symprec",
        default=DEFAULT_SYMPREC,
        type=float,
        help="Symmetry tolerance.",
    )

    return parser.parse_args()


def create_analyzer(
    structure: Structure,
    symprec: float,
):
    """
    Create a SpacegroupAnalyzer.

    Returns
    -------
    SpacegroupAnalyzer
    """

    return SpacegroupAnalyzer(
        structure,
        symprec=symprec,
    )


def analyze_structure(
    path: Path,
    symprec: float,
):
    """
    Analyze symmetry of one structure.

    Parameters
    ----------
    path

    symprec

    Returns
    -------
    dict
    """

    try:

        structure = load_structure(path)

        analyzer = create_analyzer(
            structure,
            symprec,
        )

        dataset = analyzer.get_symmetry_dataset()

        result = {
            "file": path.name,
            "space_group_number": analyzer.get_space_group_number(),
            "space_group_symbol": analyzer.get_space_group_symbol(),
            "crystal_system": analyzer.get_crystal_system(),
            "bravais_lattice": analyzer.get_lattice_type(),
            "point_group": analyzer.get_point_group_symbol(),
            "hall_symbol": dataset.hall,
            "symmetry_operations": len(
                analyzer.get_symmetry_operations()
            ),
            "symprec": symprec,
            "valid": True,
        }

        return result

    except Exception as exc:

        return {
            "file": path.name,
            "space_group_number": np.nan,
            "space_group_symbol": "Invalid",
            "crystal_system": "Invalid",
            "bravais_lattice": "Invalid",
            "point_group": "Invalid",
            "hall_symbol": "Invalid",
            "symmetry_operations": np.nan,
            "symprec": symprec,
            "valid": False,
            "error": str(exc),
        }


def print_summary(
    df: pd.DataFrame,
):
    """
    Print a dataset summary.
    """

    if not print_analysis_summary(
        df,
        title="Symmetry Analysis Summary",
        metrics=[],
        precision=3,
    ):
        return

    valid = df.loc[
        df["valid"]
    ]

    print(
        f"Unique space groups  : "
        f"{valid['space_group_symbol'].nunique()}"
    )

    print(
        f"Unique point groups  : "
        f"{valid['point_group'].nunique()}"
    )

    print(
        f"Unique Hall symbols  : "
        f"{valid['hall_symbol'].nunique()}"
    )

    print()

    print("Crystal systems")

    counts = (
        valid["crystal_system"]
        .value_counts()
        .sort_index()
    )

    total = counts.sum()

    for system, count in counts.items():

        percentage = 100 * count / total

        print(
            f"{system:15s} "
            f"{count:8d} "
            f"({percentage:5.1f}%)"
        )

    print()

    print_top_categories(
        df,
        column="space_group_symbol",
        title="Top 10 space groups",
        top_n=10,
        width=12,
    )


def save_csv(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Save analysis results.
    """

    valid = df.loc[
        df["valid"]
    ]

    extra_summary = []

    if not valid.empty:

        extra_summary.extend(
            [
                {
                    "Metric": "Unique space groups",
                    "Value": valid["space_group_symbol"].nunique(),
                },
                {
                    "Metric": "Unique point groups",
                    "Value": valid["point_group"].nunique(),
                },
                {
                    "Metric": "Unique Hall symbols",
                    "Value": valid["hall_symbol"].nunique(),
                },
                {
                    "Metric": "Unique crystal systems",
                    "Value": valid["crystal_system"].nunique(),
                },
                {
                    "Metric": "Mean symmetry operations",
                    "Value": valid["symmetry_operations"].mean(),
                },
                {
                    "Metric": "Median symmetry operations",
                    "Value": valid["symmetry_operations"].median(),
                },
                {
                    "Metric": "Maximum symmetry operations",
                    "Value": valid["symmetry_operations"].max(),
                },
            ]
        )

    save_analysis_results(
        df,
        output_dir,
        "symmetry_analysis.csv",
        [
            "space_group_number",
            "symmetry_operations",
        ],
        extra_summary=extra_summary,
    )


def create_all_plots(
    df: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    valid = df.loc[
        df["valid"]
    ]

    if valid.empty:
        return

    plot_category_distribution(
        valid["crystal_system"],
        xlabel="Count",
        ylabel="Crystal system",
        title="Crystal system distribution",
        filename="crystal_system_distribution",
        output_dir=output_dir,
        top_n=20,
    )

    plot_category_distribution(
        valid["space_group_symbol"],
        xlabel="Count",
        ylabel="Space group",
        title="Most common space groups",
        filename="space_group_distribution",
        output_dir=output_dir,
        top_n=20,
    )

    plot_histogram(
        valid["space_group_number"],
        xlabel="Space group number",
        title="Space group number distribution",
        filename="space_group_number_histogram",
        output_dir=output_dir,
        bins=230,
    )

    plot_category_distribution(
        valid["point_group"],
        xlabel="Count",
        ylabel="Point group",
        title="Point group distribution",
        filename="point_group_distribution",
        output_dir=output_dir,
        top_n=20,
    )

    plot_histogram(
        valid["symmetry_operations"],
        xlabel="Number of symmetry operations",
        title="Symmetry operations",
        filename="symmetry_operations_histogram",
        output_dir=output_dir,
        bins=30,
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
            "symmetry",
        )

        worker = partial(
            analyze_structure,
            symprec=args.symprec,
        )

        rows = parallel_map(
            compound_files,
            worker,
            args.workers,
            f"Analyzing {compound} symmetry",
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
