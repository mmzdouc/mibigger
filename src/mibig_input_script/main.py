#!/usr/bin/env python3

import argparse
from Bio import Entrez
from importlib import metadata
import json
from pathlib import Path
from typing import Dict

from mibig_input_script.aux.parse_arguments import parse_arguments
from mibig_input_script.aux.read_functions import get_curator_email
from mibig_input_script.aux.read_functions import read_mibig_json
from mibig_input_script.aux.verify_existence_entry import verify_existence_entry

from mibig_input_script.classes.mibig_entry_class import MibigEntry
from mibig_input_script.classes.class_minimal import Minimal
from mibig_input_script.classes.class_changelog import Changelog


VERSION = metadata.version("mibig-input-script")
ROOT = Path(__file__).resolve().parent
CURATION_ROUND = "next"


def get_mibig_minimal(
    existing_mibig: Dict | None, args: argparse.Namespace, ROOT: Path
) -> Dict:
    """Collect information for a minimal MIBiG entry

    Parameters:
        `existing_mibig` : existing MIBiG entry as `dict` or None
        `args` : arguments provided by user
        `ROOT` : `Path` object indicating "root" directory of script

    Returns:
        JSON-compatible dict
    """
    minimal = Minimal()

    if existing_mibig is not None:
        minimal.load_existing(existing_mibig)
    else:
        minimal.get_new_mibig_accession(args, ROOT)

    minimal.get_input()

    return minimal.export_attributes_to_dict()


def get_changelog(
    existing_mibig: Dict | None,
    args: argparse.Namespace,
    ROOT: Path,
    CURATION_ROUND: str,
) -> Dict:
    """Collect information for a minimal MIBiG entry

    Parameters:
        `existing_mibig` : existing MIBiG entry as `dict` or None
        `args` : arguments provided by user
        `ROOT` : `Path` object indicating "root" directory of script
        `CURATION_ROUND` : `str` of mibig curation round version to use in changelog

    Returns:
        JSON-compatible dict
    """
    changelog = Changelog(args.curator, CURATION_ROUND)

    print(changelog)
    # EXPAND THIS PART


def main() -> None:
    """Entry point of the program.

    Serves as main entry point for the program execution.
    It performs the following tasks:
    - Create the program interface via `argparse`
    - (Optional) Perform checks on input entry
    - Initialize a `MibigEntry()` instance
    - (Optional) Parse an existing MIBiG entry
    - Accept input to add to entry, run verification checks
    - Write the changelog
    - Dump the entry as json file

    Parameters:
        None

    Returns:
        None
    """
    args = parse_arguments(VERSION, ROOT)

    Entrez.email = get_curator_email(ROOT, args.curator)

    path_existing = verify_existence_entry(args.existing, ROOT)

    existing_mibig = read_mibig_json(path_existing)

    # add a Changelog class -> at the end -> people must confirm what they change

    # temporary
    minimal_dict = get_mibig_minimal(existing_mibig, args, ROOT)

    print(minimal_dict)

    #####

    # ~ mibig_entry = MibigEntry(get_mibig_minimal())

    # generate_minimal_entry() is the loader, which has the minimal entry class in the file (including the reading/assignment/testing parts

    # initialize the class

    # if cond to check if existing - if yes, read existing and store in object

    # get input


if __name__ == "__main__":
    main()
