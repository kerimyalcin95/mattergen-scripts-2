from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def safe_statistics(values):
    """
    Compute common statistics for a numeric sequence.

    Parameters
    ----------
    values
        Numeric sequence.

    Returns
    -------
    tuple
        mean, median, std, min, max
    """

    values = pd.Series(
        values,
        dtype="float64",
    ).dropna()

    if values.empty:
        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )

    return (
        values.mean(),
        values.median(),
        values.std(),
        values.min(),
        values.max(),
    )


def summarize_dataframe(
    df: pd.DataFrame,
    columns,
):
    """
    Compute common statistics for selected DataFrame columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.

    columns : iterable[str]
        Columns to summarize.

    Returns
    -------
    dict
        Statistics keyed by column name.
    """

    summary = {}

    for column in columns:

        values = (
            df[column]
            .dropna()
        )

        if values.empty:
            continue

        (
            mean,
            median,
            std,
            minimum,
            maximum,
        ) = safe_statistics(values)

        summary[column] = {
            "mean": mean,
            "median": median,
            "std": std,
            "min": minimum,
            "max": maximum,
        }

    return summary


def create_summary_rows(
    df: pd.DataFrame,
    columns,
):
    """
    Create summary rows suitable for CSV output.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.

    columns : iterable[str]
        Numeric columns to summarize.

    Returns
    -------
    list[dict]
        Summary rows with Metric and Value fields.
    """

    summary = [
        {
            "Metric": "Structures",
            "Value": len(df),
        },
        {
            "Metric": "Valid structures",
            "Value": int(df["valid"].sum()),
        },
        {
            "Metric": "Invalid structures",
            "Value": int((~df["valid"]).sum()),
        },
    ]

    for column in columns:

        values = (
            df[column]
            .dropna()
        )

        if values.empty:
            continue

        (
            mean,
            median,
            std,
            minimum,
            maximum,
        ) = safe_statistics(values)

        summary.extend(
            [
                {
                    "Metric": f"{column}_mean",
                    "Value": mean,
                },
                {
                    "Metric": f"{column}_median",
                    "Value": median,
                },
                {
                    "Metric": f"{column}_std",
                    "Value": std,
                },
                {
                    "Metric": f"{column}_min",
                    "Value": minimum,
                },
                {
                    "Metric": f"{column}_max",
                    "Value": maximum,
                },
            ]
        )

    return summary


def save_summary_csv(
    summary: list[dict],
    output_dir: Path,
) -> None:
    """
    Save summary rows to summary.csv.
    """

    pd.DataFrame(summary).to_csv(
        output_dir / "summary.csv",
        index=False,
    )


def save_failed_files(
    df: pd.DataFrame,
    output_dir: Path,
) -> None:
    """
    Save failed structure files to failed_files.csv.
    """

    failed = df.loc[
        df["valid"] == False
    ]

    columns = [
        "file",
    ]

    if "error" in df.columns:
        columns.append(
            "error",
        )

    failed[
        columns
    ].to_csv(
        output_dir / "failed_files.csv",
        index=False,
    )


def print_analysis_summary(
    df: pd.DataFrame,
    title: str,
    metrics,
    precision: int = 3,
) -> bool:
    """
    Print the common part of an analysis summary.

    Parameters
    ----------
    df : pandas.DataFrame
        Analysis results.

    title : str
        Summary title.

    metrics : iterable[str]
        Numeric columns to summarize.

    precision : int
        Number of decimal places.

    Returns
    -------
    bool
        True when valid structures exist.
        False when all structures are invalid.
    """

    valid = int(
        df["valid"].sum()
    )

    invalid = len(df) - valid

    print("=" * 60)
    print(title)
    print("=" * 60)
    print()

    print(
        f"Structures            : {len(df)}"
    )

    print(
        f"Successfully analyzed : {valid}"
    )

    print(
        f"Failed                : {invalid}"
    )

    print()

    if valid == 0:
        return False

    statistics = summarize_dataframe(
        df.loc[df["valid"]],
        metrics,
    )

    for column, values in statistics.items():

        print(column)

        print(
            f"  Mean   : {values['mean']:.{precision}f}"
        )

        print(
            f"  Median : {values['median']:.{precision}f}"
        )

        print(
            f"  Std    : {values['std']:.{precision}f}"
        )

        print(
            f"  Min    : {values['min']:.{precision}f}"
        )

        print(
            f"  Max    : {values['max']:.{precision}f}"
        )

        print()

    return True


def print_top_categories(
    df: pd.DataFrame,
    column: str,
    title: str,
    top_n: int = 10,
    width: int = 20,
) -> None:
    """
    Print the most common categorical values.

    Parameters
    ----------
    df : pandas.DataFrame
        Analysis results.

    column : str
        Column containing categorical values.

    title : str
        Section title.

    top_n : int
        Maximum number of categories.

    width : int
        Display width for the category.
    """

    values = df.loc[
        df["valid"],
        column,
    ].dropna()

    if values.empty:
        return

    counts = values.value_counts()

    print(title)

    for category, count in counts.head(top_n).items():

        print(
            f"{str(category):{width}s} {count}"
        )

    print()


def save_analysis_results(
    df: pd.DataFrame,
    output_dir: Path,
    filename: str,
    numeric_columns,
    extra_summary=None,
) -> None:
    """
    Save analysis results, summary statistics, and failed files.
    """

    df.to_csv(
        output_dir / filename,
        index=False,
    )

    summary = create_summary_rows(
        df,
        numeric_columns,
    )

    if extra_summary is not None:
        summary.extend(
            extra_summary,
        )

    save_summary_csv(
        summary,
        output_dir,
    )

    save_failed_files(
        df,
        output_dir,
    )


def statistics_to_dict(
    values,
    prefix: str,
    suffix: str = "",
):
    """
    Compute common statistics and return them as a dictionary.

    Parameters
    ----------
    values
        Numeric sequence.

    prefix
        Prefix used for the statistic names.

    suffix
        Optional suffix appended after the statistic name.

    Returns
    -------
    dict
        Dictionary containing mean, median, std, min and max.
    """

    (
        mean,
        median,
        std,
        minimum,
        maximum,
    ) = safe_statistics(
        values,
    )

    return {
        f"{prefix}_mean{suffix}": mean,
        f"{prefix}_median{suffix}": median,
        f"{prefix}_std{suffix}": std,
        f"{prefix}_min{suffix}": minimum,
        f"{prefix}_max{suffix}": maximum,
    }
