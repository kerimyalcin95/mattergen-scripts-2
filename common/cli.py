from __future__ import annotations

import argparse
from pathlib import Path


def create_common_parser(
    description: str,
) -> argparse.ArgumentParser:
    """
    Create an argument parser with common analysis options.

    Parameters
    ----------
    description
        Description shown by argparse.

    Returns
    -------
    argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(
        description=description,
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input folder containing crystal structures.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Base output directory.",
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search subdirectories recursively.",
    )

    parser.add_argument(
        "--workers",
        default=-1,
        type=int,
        help="-1 uses all CPU cores.",
    )

    return parser
