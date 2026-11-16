"""Shared input/output utilities for crystal structure analyses."""

from __future__ import annotations

from pathlib import Path

from ase.io import read as ase_read
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

SUPPORTED_EXTENSIONS = {
    ".cif",
    ".xyz",
    ".extxyz",
}


def find_structure_files(
    folder: Path,
    recursive: bool,
) -> list[Path]:
    """
    Find supported crystal structure files.

    Parameters
    ----------
    folder
        Input directory containing structure files.

    recursive
        Search subdirectories recursively when True.

    Returns
    -------
    list[Path]
        Sorted list of supported structure files.
    """

    iterator = (
        folder.rglob("*")
        if recursive
        else folder.glob("*")
    )

    files = [
        path
        for path in iterator
        if path.is_file()
        and path.suffix.lower() in SUPPORTED_EXTENSIONS
    ]

    files.sort()

    return files


def load_structure(
    path: Path,
) -> Structure:
    """
    Load a crystal structure using pymatgen or ASE.

    CIF files are loaded directly with pymatgen.
    XYZ and EXTXYZ files are loaded through ASE and converted
    to a pymatgen Structure.

    Parameters
    ----------
    path
        Structure file.

    Returns
    -------
    Structure
        Loaded pymatgen Structure.
    """

    if path.suffix.lower() == ".cif":
        return Structure.from_file(path)

    atoms = ase_read(path)

    return AseAtomsAdaptor.get_structure(atoms)


def ensure_output_directory(
    output_dir: Path,
) -> None:
    """
    Create an output directory if necessary.

    Parameters
    ----------
    output_dir
        Directory to create.
    """

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


def create_analysis_output_directory(
    base_output: Path,
    analysis_name: str,
) -> Path:
    """
    Create a dedicated output directory for an analysis.

    Parameters
    ----------
    base_output
        Base output directory.

    analysis_name
        Name of the analysis subdirectory.

    Returns
    -------
    Path
        Created analysis output directory.
    """

    output_dir = base_output / analysis_name

    ensure_output_directory(
        output_dir,
    )

    return output_dir
