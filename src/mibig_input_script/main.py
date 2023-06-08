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

from mibig_input_script.classes.class_write_mibig import WriteMibig
from mibig_input_script.classes.class_base import Base
from mibig_input_script.classes.class_changelog import Changelog


VERSION = metadata.version("mibig-input-script")
ROOT = Path(__file__).resolve().parent
CURATION_ROUND = "next"


def get_mibig_dict(
    existing_mibig: Dict | None, args: argparse.Namespace, ROOT: Path
) -> Dict:
    """Collect information to create a minimal MIBiG entry

    Parameters:
        `existing_mibig` : existing MIBiG entry as `dict` or None
        `args` : arguments provided by user
        `ROOT` : `Path` object indicating "root" directory of script

    Returns:
        JSON-compatible dict
    """
    mibig = Base()

    if existing_mibig is not None:
        mibig.load_existing_entry(existing_mibig)
    else:
        mibig.get_new_mibig_accession(args, ROOT)
        mibig.create_new_entry()

    mibig.get_input(existing_mibig)

    return mibig.export_attributes_to_dict()


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

    if existing_mibig is None:
        return changelog.create_new_changelog()
    elif existing_mibig["changelog"][-1]["version"] != CURATION_ROUND:
        return changelog.create_new_entry_changelog(existing_mibig)
    else:
        return changelog.append_last_entry_changelog(existing_mibig)


def write_mibig_entry(
    ROOT: Path, path_existing: Path, mibig_dict: Dict, changelog_dict: Dict
) -> None:
    """Concatenate information and create/modify a MIBiG entry

    Parameters:
        `ROOT` : `Path` object indicating "root" directory of script
        `path_existing` : `Path` object indicating location of optional existing entry
        `mibig_dict` : `dict` containing the mibig entry information
        `changelog_dict` : `dict` containing a changelog

    Returns:
        JSON-compatible dict
    """
    mibig_entry = WriteMibig()
    mibig_entry.concatenate_dicts(mibig_dict, changelog_dict)

    # CONTINUE HERE

    # call duplicate check

    # get json string for validation script (can be called from aux)

    if path_existing is None:
        new_entry_path = ROOT.joinpath("mibig_next_ver")
        mibig_entry.export_to_json(new_entry_path)
        return
    else:
        # call method to write to path
        return


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

    mibig_dict = get_mibig_dict(existing_mibig, args, ROOT)

    changelog_dict = get_changelog(existing_mibig, args, ROOT, CURATION_ROUND)

    # EXPAND HERE FOR ADDITIONAL ENTRIES

    write_mibig_entry(ROOT, path_existing, mibig_dict, changelog_dict)


if __name__ == "__main__":
    main()
