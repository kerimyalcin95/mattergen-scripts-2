#!/usr/bin/env python3
"""
Interactive launcher for MatterGen-Scripts.

The launcher provides one interface for:

- downloading reference databases,
- selecting multiple chemical systems,
- selecting existing datasets,
- selecting multiple analyses,
- configuring analysis parameters,
- configuring worker count,
- running the selected workflow.

The scientific processing and file discovery remain inside the
individual numbered scripts. This launcher only orchestrates them.
"""

from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import questionary
from questionary import Choice
from rich.console import Console
from rich.panel import Panel
from rich.table import Table


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent

DEFAULT_DATA_DIRECTORY = (
    BASE_DIR / "../../data"
).resolve()

DEFAULT_RESULTS_DIRECTORY = (
    BASE_DIR / "../../results"
).resolve()


# -----------------------------------------------------------------------------
# Console
# -----------------------------------------------------------------------------

console = Console()


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Database:
    """
    Describe a supported reference database.

    Attributes
    ----------
    name
        Human-readable database name.

    directory
        Directory name used under the data/ and results/ directories.

    identifier
        Identifier accepted by 00_download_database.py.
    """

    name: str
    directory: str
    identifier: str


DATABASES = {
    "Materials Project": Database(
        name="Materials Project",
        directory="MaterialsProject",
        identifier="mp",
    ),
    "OQMD": Database(
        name="OQMD",
        directory="OQMD",
        identifier="oqmd",
    ),
    "JARVIS": Database(
        name="JARVIS",
        directory="JARVIS",
        identifier="jarvis",
    ),
    "Alexandria": Database(
        name="Alexandria",
        directory="Alexandria",
        identifier="alexandria",
    ),
    "COD": Database(
        name="COD",
        directory="COD",
        identifier="cod",
    ),
}


MATTERGEN_DIRECTORY = "MatterGen"


ANALYSES = {
    "Basic statistics": "01_basic_statistics.py",
    "Symmetry": "02_symmetry_analysis.py",
    "Coordination": "03_coordination_analysis.py",
    "Bond": "04_bond_analysis.py",
    "RDF": "05_rdf_analysis.py",
    "Fingerprint": "06_fingerprint_analysis.py",
}


# -----------------------------------------------------------------------------
# User input
# -----------------------------------------------------------------------------


def ask_text(
    message: str,
    default: str = "",
):
    """
    Ask for text input and handle cancellation.
    """

    answer = questionary.text(
        message,
        default=default,
    ).ask()

    if answer is None:
        raise KeyboardInterrupt

    return answer.strip()


def ask_confirm(
    message: str,
    default: bool = True,
):
    """
    Ask for a yes/no confirmation and handle cancellation.
    """

    answer = questionary.confirm(
        message,
        default=default,
    ).ask()

    if answer is None:
        raise KeyboardInterrupt

    return answer


def select_action():
    """
    Select the workflow action.

    Returns
    -------
    str
        Selected action.
    """

    answer = questionary.select(
        "What do you want to do?",
        choices=[
            "Download datasets",
            "Analyze datasets",
            "Download + analyze",
            "Exit",
        ],
    ).ask()

    if answer is None:
        raise KeyboardInterrupt

    return answer


def select_databases():
    """
    Select databases to download.

    Returns
    -------
    list[Database]
        Selected database definitions.
    """

    selected = questionary.checkbox(
        "Select databases to download:",
        choices=[
            Choice(
                title=database.name,
                value=database.name,
            )
            for database in DATABASES.values()
        ],
        validate=lambda values: (
            True
            if values
            else "Select at least one database."
        ),
    ).ask()

    if selected is None:
        raise KeyboardInterrupt

    return [
        DATABASES[name]
        for name in selected
    ]


def select_analysis_datasets():
    """
    Select datasets to analyze.

    MatterGen is included as a locally generated dataset.
    """

    choices = [
        MATTERGEN_DIRECTORY,
        *DATABASES.keys(),
    ]

    selected = questionary.checkbox(
        "Select datasets to analyze:",
        choices=choices,
        validate=lambda values: (
            True
            if values
            else "Select at least one dataset."
        ),
    ).ask()

    if selected is None:
        raise KeyboardInterrupt

    return selected


def enter_chemical_systems():
    """
    Enter one or more chemical systems.

    Chemical systems are separated by whitespace.

    Examples
    --------
    Al-O Fe-O Ti-O Y-Al-O
    """

    answer = ask_text(
        "Chemical systems (e.g. Al-O Fe-O Y-Al-O):",
    )

    systems = answer.split()

    if not systems:
        raise ValueError(
            "At least one chemical system is required."
        )

    return systems


def select_analyses():
    """
    Select analyses to execute.

    Returns
    -------
    list[tuple[str, str]]
        Analysis display name and script filename.
    """

    selected = questionary.checkbox(
        "Select analyses:",
        choices=list(ANALYSES.keys()),
        validate=lambda values: (
            True
            if values
            else "Select at least one analysis."
        ),
    ).ask()

    if selected is None:
        raise KeyboardInterrupt

    return [
        (name, ANALYSES[name])
        for name in selected
    ]


def select_workers():
    """
    Select the number of parallel workers.

    -1 means all available CPU cores.
    """

    def validate_workers(value: str):
        try:
            workers = int(value)
        except ValueError:
            return "Enter an integer."

        if workers == 0:
            return "Workers cannot be zero."

        if workers < -1:
            return "Use -1 or a positive integer."

        return True

    answer = questionary.text(
        "Number of workers (-1 = all CPU cores):",
        default="-1",
        validate=validate_workers,
    ).ask()

    if answer is None:
        raise KeyboardInterrupt

    return int(answer)


def select_symprec():
    """
    Select the symmetry tolerance used by symmetry analysis.
    """

    def validate_symprec(value: str):
        try:
            number = float(value)
        except ValueError:
            return "Enter a number."

        if number <= 0:
            return "Symmetry tolerance must be greater than zero."

        return True

    answer = questionary.text(
        "Symmetry tolerance (symprec):",
        default="0.1",
        validate=validate_symprec,
    ).ask()

    if answer is None:
        raise KeyboardInterrupt

    return float(answer)


def select_rdf_parameters():
    """
    Select RDF analysis parameters.

    Returns
    -------
    tuple[float, float]
        rmax and bin_size.
    """

    def validate_positive(value: str):
        try:
            number = float(value)
        except ValueError:
            return "Enter a number."

        if number <= 0:
            return "Value must be greater than zero."

        return True

    rmax = questionary.text(
        "RDF maximum distance rmax (Å):",
        default="10.0",
        validate=validate_positive,
    ).ask()

    if rmax is None:
        raise KeyboardInterrupt

    bin_size = questionary.text(
        "RDF bin size (Å):",
        default="0.02",
        validate=validate_positive,
    ).ask()

    if bin_size is None:
        raise KeyboardInterrupt

    return (
        float(rmax),
        float(bin_size),
    )


def select_analysis_parameters(
    analyses: list[tuple[str, str]],
):
    """
    Select parameters required by the selected analysis scripts.

    Returns
    -------
    dict
        Script-specific analysis parameters.
    """

    parameters = {}

    selected_names = {
        name
        for name, _ in analyses
    }

    if "Symmetry" in selected_names:

        parameters["symprec"] = (
            select_symprec()
        )

    if "RDF" in selected_names:

        (
            rmax,
            bin_size,
        ) = select_rdf_parameters()

        parameters["rmax"] = rmax
        parameters["bin_size"] = bin_size

    return parameters


def select_directories():
    """
    Select the data and results root directories.

    Returns
    -------
    tuple[Path, Path]
        Data and results directories.
    """

    data_directory = questionary.path(
        "Data directory:",
        default=str(DEFAULT_DATA_DIRECTORY),
    ).ask()

    if data_directory is None:
        raise KeyboardInterrupt

    results_directory = questionary.path(
        "Results directory:",
        default=str(DEFAULT_RESULTS_DIRECTORY),
    ).ask()

    if results_directory is None:
        raise KeyboardInterrupt

    return (
        Path(data_directory).expanduser().resolve(),
        Path(results_directory).expanduser().resolve(),
    )


# -----------------------------------------------------------------------------
# Dataset paths
# -----------------------------------------------------------------------------


def get_dataset_input_directory(
    dataset: str,
    chemical_system: str,
    data_directory: Path,
):
    """
    Return the input directory for a dataset.

    MatterGen stores generated CIF structures in:

        data/MatterGen/<chemical-system>/generated_crystals_cif

    Reference databases store structures in:

        data/<database>/<chemical-system>
    """

    if dataset == MATTERGEN_DIRECTORY:

        return (
            data_directory
            / MATTERGEN_DIRECTORY
            / chemical_system
            / "generated_crystals_cif"
        )

    database = DATABASES[dataset]

    return (
        data_directory
        / database.directory
        / chemical_system
    )


def get_dataset_output_directory(
    dataset: str,
    chemical_system: str,
    results_directory: Path,
):
    """
    Return the base results directory for a dataset.
    """

    if dataset == MATTERGEN_DIRECTORY:
        directory = MATTERGEN_DIRECTORY

    else:
        directory = DATABASES[dataset].directory

    return (
        results_directory
        / directory
        / chemical_system
    )


# -----------------------------------------------------------------------------
# Command construction
# -----------------------------------------------------------------------------


def build_download_command(
    database: Database,
    chemical_systems: list[str],
    data_directory: Path,
):
    """
    Build the command for 00_download_database.py.

    The downloader expects identifiers such as:

        mp
        oqmd
        jarvis
        alexandria
        cod
    """

    return [
        sys.executable,
        str(BASE_DIR / "00_download_database.py"),
        "--database",
        database.identifier,
        "--chemsys",
        *chemical_systems,
        "--output",
        str(data_directory),
    ]


def build_analysis_command(
    script_name: str,
    input_directory: Path,
    output_directory: Path,
    workers: int,
    analysis_parameters: dict,
):
    """
    Build the complete command for an analysis script.

    Every analysis receives the common arguments:

        --input
        --output
        --workers

    Script-specific arguments are added when required:

        02_symmetry_analysis.py
            --symprec

        05_rdf_analysis.py
            --rmax
            --bin-size

    Recursive discovery is intentionally not passed. Each analysis
    operates on the files directly inside its supplied input directory.
    """

    command = [
        sys.executable,
        str(BASE_DIR / script_name),
        "--input",
        str(input_directory),
        "--output",
        str(output_directory),
        "--workers",
        str(workers),
    ]

    if script_name == "02_symmetry_analysis.py":

        command.extend(
            [
                "--symprec",
                str(
                    analysis_parameters["symprec"]
                ),
            ]
        )

    elif script_name == "05_rdf_analysis.py":

        command.extend(
            [
                "--rmax",
                str(
                    analysis_parameters["rmax"]
                ),
                "--bin-size",
                str(
                    analysis_parameters["bin_size"]
                ),
            ]
        )

    return command


# -----------------------------------------------------------------------------
# Command execution
# -----------------------------------------------------------------------------


def format_command(
    command: list[str],
):
    """
    Format a subprocess command for terminal display.
    """

    formatted = []

    for argument in command:

        if any(
            character.isspace()
            for character in argument
        ):
            formatted.append(
                f'"{argument}"'
            )

        else:
            formatted.append(
                argument
            )

    return " ".join(
        formatted
    )


def run_command(
    command: list[str],
    title: str,
):
    """
    Execute one external command.

    Returns
    -------
    bool
        True if the command exits successfully.
    """

    console.print()

    console.print(
        Panel(
            f"[bold]{title}[/bold]",
            expand=False,
        )
    )

    console.print(
        f"[dim]$ {format_command(command)}[/dim]"
    )

    console.print()

    try:

        result = subprocess.run(
            command,
            cwd=BASE_DIR,
            check=False,
        )

    except KeyboardInterrupt:

        console.print(
            "\n[yellow]Operation cancelled.[/yellow]"
        )

        return False

    if result.returncode == 0:

        console.print(
            f"[green]✓ {title} completed[/green]"
        )

        return True

    console.print(
        f"[red]✗ {title} failed "
        f"(exit code {result.returncode})[/red]"
    )

    return False


# -----------------------------------------------------------------------------
# Workflow
# -----------------------------------------------------------------------------


def run_downloads(
    databases: list[Database],
    chemical_systems: list[str],
    data_directory: Path,
):
    """
    Download all selected databases.

    Each database receives all selected chemical systems in one
    invocation.
    """

    success = True

    for database in databases:

        command = build_download_command(
            database=database,
            chemical_systems=chemical_systems,
            data_directory=data_directory,
        )

        completed = run_command(
            command=command,
            title=f"Downloading {database.name}",
        )

        if completed:
            continue

        success = False

        continue_running = ask_confirm(
            "Download failed. Continue with the remaining databases?",
            default=True,
        )

        if not continue_running:
            return False

    return success


def run_analyses(
    datasets: list[str],
    chemical_systems: list[str],
    analyses: list[tuple[str, str]],
    data_directory: Path,
    results_directory: Path,
    workers: int,
    analysis_parameters: dict,
):
    """
    Run all selected analyses for all selected dataset/system pairs.

    The launcher does not inspect the input directories and does not
    search for structure files. The individual analysis scripts are
    responsible for file discovery and empty-directory handling.
    """

    success = True

    for dataset in datasets:

        for chemical_system in chemical_systems:

            input_directory = (
                get_dataset_input_directory(
                    dataset=dataset,
                    chemical_system=chemical_system,
                    data_directory=data_directory,
                )
            )

            output_directory = (
                get_dataset_output_directory(
                    dataset=dataset,
                    chemical_system=chemical_system,
                    results_directory=results_directory,
                )
            )

            for analysis_name, script_name in analyses:

                script_path = (
                    BASE_DIR / script_name
                )

                if not script_path.is_file():

                    console.print(
                        "[red]Missing analysis script:[/red] "
                        f"{script_path}"
                    )

                    success = False

                    continue

                command = build_analysis_command(
                    script_name=script_name,
                    input_directory=input_directory,
                    output_directory=output_directory,
                    workers=workers,
                    analysis_parameters=analysis_parameters,
                )

                completed = run_command(
                    command=command,
                    title=(
                        f"{analysis_name} | "
                        f"{dataset} | "
                        f"{chemical_system}"
                    ),
                )

                if not completed:
                    success = False

    return success


# -----------------------------------------------------------------------------
# Configuration display
# -----------------------------------------------------------------------------


def print_configuration(
    action: str,
    databases: list[Database],
    datasets: list[str],
    chemical_systems: list[str],
    analyses: list[tuple[str, str]],
    workers: int,
    analysis_parameters: dict,
    data_directory: Path,
    results_directory: Path,
):
    """
    Display the complete workflow configuration.
    """

    table = Table(
        title="Configuration",
        show_header=True,
    )

    table.add_column(
        "Setting",
        style="bold",
    )

    table.add_column(
        "Value",
    )

    database_text = ", ".join(
        database.name
        for database in databases
    )

    dataset_text = ", ".join(
        datasets
    )

    analysis_text = ", ".join(
        name
        for name, _ in analyses
    )

    chemical_system_text = ", ".join(
        chemical_systems
    )

    table.add_row(
        "Action",
        action,
    )

    table.add_row(
        "Download databases",
        database_text or "None",
    )

    table.add_row(
        "Analysis datasets",
        dataset_text or "None",
    )

    table.add_row(
        "Chemical systems",
        chemical_system_text,
    )

    table.add_row(
        "Analyses",
        analysis_text or "None",
    )

    table.add_row(
        "Workers",
        str(workers),
    )

    if "symprec" in analysis_parameters:

        table.add_row(
            "Symmetry tolerance",
            str(
                analysis_parameters["symprec"]
            ),
        )

    if "rmax" in analysis_parameters:

        table.add_row(
            "RDF rmax",
            str(
                analysis_parameters["rmax"]
            ),
        )

    if "bin_size" in analysis_parameters:

        table.add_row(
            "RDF bin size",
            str(
                analysis_parameters["bin_size"]
            ),
        )

    table.add_row(
        "Data directory",
        str(data_directory),
    )

    table.add_row(
        "Results directory",
        str(results_directory),
    )

    console.print()
    console.print(table)
    console.print()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    """
    Run the interactive MatterGen-Scripts launcher.
    """

    console.clear()

    console.print(
        Panel(
            "[bold]MatterGen-Scripts[/bold]\n"
            "Dataset download and analysis launcher",
            expand=False,
        )
    )

    console.print()

    action = select_action()

    if action == "Exit":

        console.print(
            "[dim]Goodbye.[/dim]"
        )

        return

    console.print()

    databases = []

    if action in {
        "Download datasets",
        "Download + analyze",
    }:

        databases = select_databases()

    console.print()

    chemical_systems = (
        enter_chemical_systems()
    )

    console.print()

    datasets = []
    analyses = []

    if action in {
        "Analyze datasets",
        "Download + analyze",
    }:

        datasets = (
            select_analysis_datasets()
        )

        console.print()

        analyses = select_analyses()

    console.print()

    workers = select_workers()

    console.print()

    analysis_parameters = (
        select_analysis_parameters(
            analyses
        )
    )

    console.print()

    (
        data_directory,
        results_directory,
    ) = select_directories()

    console.print()

    print_configuration(
        action=action,
        databases=databases,
        datasets=datasets,
        chemical_systems=chemical_systems,
        analyses=analyses,
        workers=workers,
        analysis_parameters=analysis_parameters,
        data_directory=data_directory,
        results_directory=results_directory,
    )

    if not ask_confirm(
        "Start workflow?",
        default=True,
    ):

        console.print(
            "[yellow]Cancelled.[/yellow]"
        )

        return

    if action in {
        "Download datasets",
        "Download + analyze",
    }:

        data_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    if action in {
        "Analyze datasets",
        "Download + analyze",
    }:

        results_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    success = True

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    if action in {
        "Download datasets",
        "Download + analyze",
    }:

        download_success = run_downloads(
            databases=databases,
            chemical_systems=chemical_systems,
            data_directory=data_directory,
        )

        success = (
            success
            and download_success
        )

    # -------------------------------------------------------------------------
    # Analysis
    # -------------------------------------------------------------------------

    if action in {
        "Analyze datasets",
        "Download + analyze",
    }:

        analysis_success = run_analyses(
            datasets=datasets,
            chemical_systems=chemical_systems,
            analyses=analyses,
            data_directory=data_directory,
            results_directory=results_directory,
            workers=workers,
            analysis_parameters=analysis_parameters,
        )

        success = (
            success
            and analysis_success
        )

    console.print()

    if success:

        console.print(
            Panel(
                "[bold green]"
                "Workflow completed successfully."
                "[/bold green]",
                expand=False,
            )
        )

    else:

        console.print(
            Panel(
                "[bold yellow]"
                "Workflow completed with errors or skipped items."
                "[/bold yellow]",
                expand=False,
            )
        )


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        console.print(
            "\n[yellow]Launcher cancelled.[/yellow]"
        )

    except Exception as error:

        console.print(
            Panel(
                f"[bold red]Launcher error[/bold red]\n\n"
                f"{type(error).__name__}: {error}",
                expand=False,
            )
        )

        raise