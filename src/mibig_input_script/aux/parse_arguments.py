"""
Parse command line arguments.

Parses command line arguments for `mibig-input-script.py`.
"""

import argparse
from pathlib import Path
from typing import List

from mibig_input_script.aux.read_functions import get_curators


def add_parser_args(VERSION: str, ROOT: Path) -> argparse.ArgumentParser:
    """Add arguments to argparse object.

    Parameters:
        VERSION: str indicating the version
        `ROOT`: `Path` object indicating "root" directory of script

    Returns:
        `argparse.ArgumentParser` object
    """
    curators = get_curators(ROOT)

    descr = (
        f"mibig-input-script v{VERSION}: " "Creation and manipulation of MIBiG entries"
    )
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument(
        "curator",
        choices=curators,
        help="Give the initials of the curator creating/manipulating the MIBiG entry.",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-n", "--new", action="store_true", help="Create a new MIBiG entry."
    )
    group.add_argument(
        "-e",
        "--existing",
        nargs=1,
        type=str,
        help="Manipulate an existing entry by specifying its MIBiG accession number.",
    )
    return parser


def parse_arguments(VERSION: str, ROOT: Path) -> argparse.Namespace:
    """Create program interface using argparse.

    Args:
        `VERSION`: str indicating the version
        `ROOT`: `Path` object indicating "root" directory of script

    Returns:
        `args` : Returns argparse object for further use
    """
    parser = add_parser_args(VERSION, ROOT)
    args = parser.parse_args()

    return args
