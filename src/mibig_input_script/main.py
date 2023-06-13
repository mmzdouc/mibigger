#!/usr/bin/env python3

import argparse
from Bio import Entrez
from importlib import metadata
import json
import jsonschema
from pathlib import Path
from typing import Dict

from mibig_input_script.aux.parse_arguments import parse_arguments
from mibig_input_script.aux.read_functions import get_curator_email
from mibig_input_script.aux.read_functions import read_mibig_json
from mibig_input_script.aux.verify_existence_entry import verify_existence_entry

from mibig_input_script.classes.class_mibig_entry import MibigEntry
from mibig_input_script.classes.class_changelog import Changelog
from mibig_input_script.classes.class_write_mibig import WriteMibig


VERSION = metadata.version("mibig-input-script")
ROOT = Path(__file__).resolve().parent
CURATION_ROUND = "next"


def get_mibig_entry(
    existing_mibig: Dict | None,
    args: argparse.Namespace,
    ROOT: Path,
    CURATION_ROUND: str,
) -> Dict:
    """Collect information to create a minimal MIBiG entry

    Parameters:
        `existing_mibig` : existing MIBiG entry as `dict` or None
        `args` : arguments provided by user
        `ROOT` : `Path` object indicating "root" directory of script
        `CURATION_ROUND` : `str` of mibig curation round version to use in changelog

    Returns:
        new/modified mibig entry as `dict`
    """
    mibig_entry = MibigEntry()

    if existing_mibig is not None:
        mibig_entry.load_existing_entry(existing_mibig)
    else:
        mibig_entry.get_new_mibig_accession(args, ROOT)
        mibig_entry.create_new_entry(args.curator, CURATION_ROUND)

    mibig_entry.get_input()

    return mibig_entry.export_attributes_to_dict()


def get_changelog(
    mibig_entry: Dict,
    args: argparse.Namespace,
    CURATION_ROUND: str,
) -> Dict:
    """Collect information for a minimal MIBiG entry

    Parameters:
        `mibig_entry` : MIBiG entry as `dict`
        `args` : arguments provided by user
        `CURATION_ROUND` : `str` of mibig curation round version to use in changelog

    Returns:
        new/modified mibig entry as `dict`
    """
    modify_changelog = Changelog(args.curator, CURATION_ROUND)

    if mibig_entry["changelog"][-1]["version"] != CURATION_ROUND:
        return modify_changelog.create_new_entry_changelog(mibig_entry)
    else:
        return modify_changelog.append_last_entry_changelog(mibig_entry)


def export_mibig_entry(mibig_entry: Dict, path_existing: Path, ROOT: Path) -> None:
    """Concatenate information and create/modify a MIBiG entry

    Parameters:
        `mibig_dict` : `dict` containing the mibig entry information
        `path_existing` : `Path` object indicating location of optional existing entry
        `ROOT` : `Path` object indicating "root" directory of script

    Returns
        None
    """
    write_mibig = WriteMibig(mibig_entry)

    write_mibig.test_duplicate_entries(ROOT)

    json_string = write_mibig.return_json_string()

    with open(ROOT.joinpath("schema.json")) as schema_handle:
        schema = json.load(schema_handle)
    try:
        jsonschema.validate(json.loads(json_string), schema)
    except jsonschema.ValidationError as e:
        print("MIBiG JSON is invalid. Error:", e)
        print("Abort file storage.")
        exit()

    if path_existing is None:
        new_entry_path = (
            ROOT.joinpath("mibig_next_ver")
            .joinpath(write_mibig.export_dict["cluster"]["mibig_accession"])
            .with_suffix(".json")
        )
        write_mibig.export_to_json(new_entry_path)
        print("New entry successfully created.")
        write_mibig.append_to_csv_existing(ROOT)
    else:
        write_mibig.export_to_json(path_existing)
        print("Existing entry successfully modified.")

    return


def main() -> None:
    """Entry point of the program.

    Serves as main entry point for the program execution.
    It performs the following tasks:
    - Create the program interface via `argparse`
    - Run auxilliary functions
    - Initialize a `MibigEntry` instance, take and test input data
    - Initialized a `Changelog` instance, write changelog entry
    - Initialize a `WriteMibig` instance, prepare for export
    - Validate resulting JSON file
    - Store JSON file

    Parameters:
        None

    Returns:
        None
    """
    args = parse_arguments(VERSION, ROOT)

    Entrez.email = get_curator_email(ROOT, args.curator)
    path_existing = verify_existence_entry(args.existing, ROOT)
    existing_mibig = read_mibig_json(path_existing)

    mibig_entry = get_mibig_entry(existing_mibig, args, ROOT, CURATION_ROUND)

    # EXPAND HERE FOR ADDITIONAL DATA ENTRIES

    mibig_entry = get_changelog(mibig_entry, args, CURATION_ROUND)

    export_mibig_entry(mibig_entry, path_existing, ROOT)


if __name__ == "__main__":
    main()
