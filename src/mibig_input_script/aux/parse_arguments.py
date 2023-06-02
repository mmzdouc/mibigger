#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import List


def get_curators() -> List[str]:
    """Get the list of possible curators defined in ../curators.csv

    Returns:
        Returns a list of curator ids
    """
    parent_dir = Path(__file__).resolve().parent.parent
    csv_path = parent_dir / ".curators.csv"

    with open(csv_path, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        return [row[0] for row in csvreader]


def add_parser_args(VERSION: str) -> argparse.Namespace:
    """Add arguments to argparse object

    Parameters:
        VERSION: str indicating the version

    Returns:
        argparse object

    Notes:
        Additional command line parameters can be easily specified
            by adding arguments below.

        Technical:
        `choices` allows to restrict options
        `action=store_true` returns True if flag is set
        `nargs` limits n of args to max 1
    """
    curators = get_curators()

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

    return parser


def parse_arguments(VERSION: str) -> argparse.Namespace:
    """Create program interface using argparse

    Args:
        VERSION: str indicating the version

    Returns:
        Returns argparse object for further use

    Notes:
        At least one of the optional flags
    """

    parser = add_parser_args(VERSION)

    return parser.parse_args()
