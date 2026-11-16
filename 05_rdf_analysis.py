#!/usr/bin/env python3
"""
05_rdf_analysis.py

Analyze radial distribution functions (RDF) of crystal structures.

Features
--------
- Reads CIF, XYZ and EXTXYZ files
- Computes radial distribution functions
- Computes:
    * First RDF peak position
    * First RDF peak height
    * RDF integral
    * Maximum RDF
- Saves:
    * CSV with one row per structure
    * Average RDF CSV
    * Summary CSV
    * Failed files CSV
    * Publication-quality PDF and PNG figures
- Prints a dataset summary

Supported libraries
-------------------
- pymatgen
- ase
- numpy
- pandas
- scipy
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

from scipy.integrate import trapezoid
from scipy.signal import find_peaks

from pymatgen.core import Structure

from pymatviz.rdf.helpers import calculate_rdf

from functools import partial

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
)

from common.cli import create_common_parser

from common.plotting import (
    plot_histogram,
    save_figure,
    configure_plotting
)

DEFAULT_R_MAX = 10.0
DEFAULT_BIN_SIZE = 0.02

configure_plotting()


def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = create_common_parser(
        "Analyze radial distribution functions."
    )

    parser.add_argument(
        "--rmax",
        default=DEFAULT_R_MAX,
        type=float,
        help="Maximum RDF distance in Å.",
    )

    parser.add_argument(
        "--bin-size",
        default=DEFAULT_BIN_SIZE,
        type=float,
        help="Histogram bin width in Å.",
    )

    return parser.parse_args()


def compute_rdf(
    structure: Structure,
    rmax: float,
    bin_size: float,
):
    """
    Compute the radial distribution function using pymatviz.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
    """

    n_bins = max(
        1,
        int(rmax / bin_size),
    )

    r, rdf = calculate_rdf(
        structure=structure,
        cutoff=rmax,
        n_bins=n_bins,
        pbc=(1, 1, 1),
    )

    return (
        np.asarray(r, dtype=float),
        np.asarray(rdf, dtype=float),
    )


def analyze_structure(
    path: Path,
    rmax: float,
    bin_size: float,
):
    """
    Analyze one crystal structure.

    Parameters
    ----------
    path

    rmax

    bin_size

    Returns
    -------
    tuple
        analysis dictionary
        RDF dataframe
    """

    try:

        structure = load_structure(
            path,
        )

        r, rdf = compute_rdf(
            structure,
            rmax,
            bin_size,
        )

        statistics = extract_rdf_statistics(
            r,
            rdf,
        )

        rdf_frame = pd.DataFrame(
            {
                "file": path.name,
                "radius": r,
                "rdf": rdf,
            }
        )

        result = create_structure_result(
            path,
            structure,
            statistics,
        )

        return (
            result,
            rdf_frame,
        )

    except Exception as exc:

        result = create_failed_result(
            path,
            {
                "first_peak_position": np.nan,
                "first_peak_height": np.nan,
                "coordination_shell_radius": np.nan,
                "rdf_integral": np.nan,
                "rdf_maximum": np.nan,
            },
            exc,
        )

        return (
            result,
            None,
        )


def print_summary(
    df: pd.DataFrame,
):
    """
    Print dataset summary.
    """

    if not print_analysis_summary(
        df,
        title="RDF Analysis Summary",
        metrics=[
            "first_peak_position",
            "first_peak_height",
            "rdf_integral",
            "rdf_maximum",
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


def extract_rdf_statistics(
    r,
    rdf,
):
    """
    Extract numerical RDF descriptors.
    """

    if len(rdf) == 0:
        return {
            "first_peak_position": np.nan,
            "first_peak_height": np.nan,
            "coordination_shell_radius": np.nan,
            "rdf_integral": np.nan,
            "rdf_maximum": np.nan,
        }

    peaks, _ = find_peaks(
        rdf,
        prominence=0.05,
    )

    if len(peaks):
        peak = peaks[0]

        peak_position = float(r[peak])
        peak_height = float(rdf[peak])

    else:
        peak_position = np.nan
        peak_height = np.nan

    return {
        "first_peak_position": peak_position,
        "first_peak_height": peak_height,
        "coordination_shell_radius": peak_position,
        "rdf_integral": float(trapezoid(rdf, r)),
        "rdf_maximum": float(np.max(rdf)),
    }


def save_csv(
    df: pd.DataFrame,
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Save CSV outputs.
    """

    save_analysis_results(
        df,
        output_dir,
        "rdf_analysis.csv",
        [
            "first_peak_position",
            "first_peak_height",
            "coordination_shell_radius",
            "rdf_integral",
            "rdf_maximum",
        ],
    )

    rdf_curves.to_csv(
        output_dir / "rdf_curves.csv",
        index=False,
    )


def plot_average_rdf(
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot the average RDF over all structures.
    """

    if rdf_curves.empty:
        return

    average = (
        rdf_curves
        .groupby("radius")["rdf"]
        .mean()
        .reset_index()
    )

    fig, ax = plt.subplots()

    ax.plot(
        average["radius"],
        average["rdf"],
        linewidth=2,
    )

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("g(r)")
    ax.set_title("Average Radial Distribution Function")

    save_figure(
        fig,
        output_dir,
        "average_rdf",
    )


def plot_rdf_heatmap(
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Plot RDF heatmap.
    """

    if rdf_curves.empty:
        return

    pivot = rdf_curves.pivot(
        index="file",
        columns="radius",
        values="rdf",
    )

    fig, ax = plt.subplots(
        figsize=(10, 8),
    )

    sns.heatmap(
        pivot,
        cmap="viridis",
        ax=ax,
        cbar_kws={
            "label": "g(r)"
        },
    )

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Structure")
    ax.set_title("Radial Distribution Functions")

    save_figure(
        fig,
        output_dir,
        "rdf_heatmap",
    )


def create_all_plots(
    df: pd.DataFrame,
    rdf_curves: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate all figures.
    """

    valid = df[
        df["valid"]
    ]

    if valid.empty:
        return

    plot_average_rdf(
        rdf_curves,
        output_dir,
    )

    plot_rdf_heatmap(
        rdf_curves,
        output_dir,
    )

    plot_histogram(
        valid["first_peak_position"],
        xlabel="Distance (Å)",
        title="First RDF Peak Position",
        filename="first_peak_distribution",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["first_peak_height"],
        xlabel="g(r)",
        title="First RDF Peak Height",
        filename="first_peak_height_distribution",
        output_dir=output_dir,
    )

    plot_histogram(
        valid["coordination_shell_radius"],
        xlabel="Distance (Å)",
        title="First Coordination Shell",
        filename="coordination_shell_distribution",
        output_dir=output_dir,
    )


def main():

    args = parse_arguments()

    output_dir = create_analysis_output_directory(
        args.output,
        "rdf",
    )

    print()
    print("Searching for structures...")

    files = find_structure_files(
        args.input,
        args.recursive,
    )

    print(
        f"Found {len(files)} structures."
    )

    if len(files) == 0:

        print(
            "No supported files found."
        )

        return

    print()

    worker = partial(
        analyze_structure,
        rmax=args.rmax,
        bin_size=args.bin_size,
    )

    results = parallel_map(
        files,
        worker,
        args.workers,
        "Computing RDF",
    )

    analysis_rows = []
    rdf_frames = []

    for row, rdf in results:

        analysis_rows.append(
            row,
        )

        if rdf is not None:
            rdf_frames.append(
                rdf,
            )

    df = pd.DataFrame(
        analysis_rows,
    )

    if rdf_frames:

        rdf_curves = pd.concat(
            rdf_frames,
            ignore_index=True,
        )

    else:

        rdf_curves = pd.DataFrame(
            columns=[
                "file",
                "radius",
                "rdf",
            ]
        )

    print_summary(
        df,
    )

    save_csv(
        df,
        rdf_curves,
        output_dir,
    )

    print(
        "Generating plots..."
    )

    create_all_plots(
        df,
        rdf_curves,
        output_dir,
    )

    print()

    print(
        f"Results written to {output_dir.resolve()}"
    )

    print(
        "Done."
    )


if __name__ == "__main__":
    main()
