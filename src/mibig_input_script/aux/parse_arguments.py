#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List

from mibig_input_script.aux.read_functions import get_curators


def add_parser_args(VERSION: str, ROOT: Path) -> argparse.ArgumentParser:
    """Add arguments to argparse object

    Parameters:
        VERSION: str indicating the version
        `ROOT`: `Path` object indicating "root" directory of script

    Returns:
        `argparse.ArgumentParser` object

    Notes:
        Additional command line parameters can be easily specified
            by adding arguments below.

        Technical:
        `choices` allows to restrict options
        `action=store_true` returns True if flag is set
        `nargs` limits n of args to max 1
    """
    curators = get_curators(ROOT)

    parser = argparse.ArgumentParser(
        description=f"""
            mibig-input-script v{VERSION}:
            Creation and manipulation of MIBiG entries
        """
    )

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

    parser.add_argument(
        "--minimal",
        action="store_true",
        help="Create/manipulate the minimal information of a MIBiG entry.",
    )

    parser.add_argument(
        "--gene_ann",
        action="store_true",
        help="TBA Add information on domains, functions, products, ... in locus.",
    )

    parser.add_argument(
        "--ripp_ann",
        action="store_true",
        help="TBA Add information on leader/core peptides, crosslinks, ... of RiPPs.",
    )

    parser.add_argument(
        "--all",
        action="store_true",
        help="TBA Add full information suite at once",
    )

    return parser


def parse_arguments(VERSION: str, ROOT: Path) -> argparse.Namespace:
    """Create program interface using argparse

    Args:
        `VERSION`: str indicating the version
        `ROOT`: `Path` object indicating "root" directory of script

    Returns:
        `args` : Returns argparse object for further use
    """

    parser = add_parser_args(VERSION, ROOT)
    args = parser.parse_args()

    return args
