from __future__ import annotations

from pathlib import Path

import numpy as np

from pymatgen.core import Structure


def create_structure_result(
    path: Path,
    structure: Structure,
    values: dict,
) -> dict:
    """
    Create a successful structure-analysis result.

    Parameters
    ----------
    path
        Structure file.

    structure
        Loaded pymatgen structure.

    values
        Analysis-specific result values.

    Returns
    -------
    dict
        Result containing common structure metadata and
        analysis-specific values.
    """

    result = {
        "file": path.name,
        "formula": structure.composition.formula,
        "reduced_formula": structure.composition.reduced_formula,
        "num_atoms": len(structure),
    }

    result.update(values)

    result["valid"] = True

    return result


def create_failed_result(
    path: Path,
    values: dict,
    error: Exception,
) -> dict:
    """
    Create a failed structure-analysis result.

    Parameters
    ----------
    path
        Structure file.

    values
        Analysis-specific fields and their failure values.

    error
        Exception raised during analysis.

    Returns
    -------
    dict
        Failed result containing common fields and
        analysis-specific failure values.
    """

    result = {
        "file": path.name,
        "formula": None,
        "reduced_formula": None,
        "num_atoms": np.nan,
    }

    result.update(values)

    result["valid"] = False
    result["error"] = str(error)

    return result