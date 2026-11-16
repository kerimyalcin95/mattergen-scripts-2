from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def configure_plotting() -> None:
    """
    Configure the project-wide Matplotlib and Seaborn settings.
    """

    sns.set_theme(
        style="whitegrid",
    )

    plt.rcParams["figure.figsize"] = (8, 5)
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["axes.titlesize"] = 14
    plt.rcParams["axes.labelsize"] = 12


def save_figure(
    figure,
    output_dir: Path,
    filename: str,
) -> None:
    """
    Save a Matplotlib figure as PNG and PDF.
    """

    figure.tight_layout()

    figure.savefig(
        output_dir / f"{filename}.png",
        dpi=300,
    )

    figure.savefig(
        output_dir / f"{filename}.pdf",
    )

    plt.close(
        figure,
    )


def plot_histogram(
    values,
    xlabel: str,
    title: str,
    filename: str,
    output_dir: Path,
    bins: int = 40,
) -> None:
    """
    Plot a histogram and save it as PNG and PDF.
    """

    values = pd.Series(
        values,
    ).dropna()

    if values.empty:
        return

    figure, axis = plt.subplots()

    sns.histplot(
        values,
        bins=bins,
        kde=True,
        stat="count",
        ax=axis,
    )

    axis.set_xlabel(
        xlabel,
    )

    axis.set_ylabel(
        "Count",
    )

    axis.set_title(
        title,
    )

    save_figure(
        figure,
        output_dir,
        filename,
    )


def plot_category_distribution(
    values,
    xlabel: str,
    ylabel: str,
    title: str,
    filename: str,
    output_dir: Path,
    top_n: int = 20,
) -> None:
    """
    Plot the most common categorical values.

    Parameters
    ----------
    values
        Categorical values to count.

    xlabel
        Label for the x-axis.

    ylabel
        Label for the y-axis.

    title
        Figure title.

    filename
        Output filename without extension.

    output_dir
        Directory where the figure is saved.

    top_n
        Maximum number of categories to display.
    """

    values = pd.Series(values).dropna()

    if values.empty:
        return

    counts = (
        values
        .value_counts()
        .head(top_n)
    )

    if counts.empty:
        return

    figure, axis = plt.subplots(
        figsize=(10, 6),
    )

    sns.barplot(
        x=counts.values,
        y=counts.index,
        ax=axis,
    )

    axis.set_xlabel(
        xlabel,
    )

    axis.set_ylabel(
        ylabel,
    )

    axis.set_title(
        title,
    )

    save_figure(
        figure,
        output_dir,
        filename,
    )
