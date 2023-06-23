#!/usr/bin/env python3

import argparse
from Bio import Entrez
from importlib import metadata
from pathlib import Path
import readline
from typing import Dict

from mibig_input_script.aux.parse_arguments import parse_arguments
from mibig_input_script.aux.read_functions import get_curator_email
from mibig_input_script.aux.read_functions import read_mibig_json
from mibig_input_script.aux.verify_existence_entry import verify_existence_entry
from mibig_input_script.aux.validation_mibig import validation_mibig
from mibig_input_script.aux.read_functions import get_curators

from mibig_input_script.classes.class_ripp import Ripp
from mibig_input_script.classes.class_mibig_entry import MibigEntry
from mibig_input_script.classes.class_changelog import Changelog
from mibig_input_script.classes.class_write_mibig import WriteMibig
from mibig_input_script.classes.class_genes import Genes


VERSION = metadata.version("mibig-input-script")
ROOT = Path(__file__).resolve().parent
CURATION_ROUND = "next"


def get_mibig_entry(
    existing_mibig: Dict | None,
    args: argparse.Namespace,
    ROOT: Path,
    CURATION_ROUND: str,
) -> Dict:
    """Collect information to create a minimal MIBiG entry.

    Parameters:
        `existing_mibig` : existing MIBiG entry as `dict` or None
        `args` : arguments provided by user
        `ROOT` : `Path` object indicating "root" directory of script
        `CURATION_ROUND` : `str` of mibig curation round version to use in changelog

    Returns:
        new/modified mibig entry as `dict`
    """
    mibig_entry_object = MibigEntry()

    if existing_mibig is not None:
        mibig_entry_object.load_existing_entry(existing_mibig)
    else:
        curator_id = get_curators(ROOT)[args.curator]
        mibig_entry_object.get_new_mibig_accession(curator_id, ROOT)
        mibig_entry_object.create_new_entry(curator_id, CURATION_ROUND)

    mibig_entry_object.get_input()

    return mibig_entry_object.export_dict()


def get_ripp(mibig_entry: Dict, ROOT: Path) -> Dict:
    """Collect data related to RiPPs.

    Parameters:
        `mibig_entry` : MIBiG entry as `dict`
        `ROOT` : A Path object pointing towards the root dir

    Returns:
        new/modified mibig entry as `dict`
    """
    if "RiPP" not in mibig_entry["cluster"]["biosyn_class"]:
        print("Cannot add information if 'biosynthetic class' is not RiPP!")
        return mibig_entry
    else:
        ripp_object = Ripp(mibig_entry, ROOT)

        try:
            ripp_object.mibig_dict["cluster"]["ripp"]
        except KeyError:
            ripp_object.initialize_ripp_entry()

        ripp_object.get_input()

        return ripp_object.export_dict()


def get_genes(mibig_entry: Dict, ROOT: Path) -> Dict:
    """Collect data related to gene.

    Parameters:
        `mibig_entry` : MIBiG entry as `dict`
        `ROOT` : A Path object pointing towards the "root" dir

    Returns:
        new/modified mibig entry as `dict`
    """
    genes_object = Genes(mibig_entry, ROOT)

    try:
        genes_object.mibig_dict["cluster"]["genes"]
    except KeyError:
        genes_object.initialize_genes_entry()

    genes_object.get_gene_info()

    return genes_object.export_dict()


def get_optional_data(mibig_entry: Dict, ROOT: Path) -> Dict:
    """Menu to add additional (optional) data.

    Parameters:
        `mibig_entry` : MIBiG entry as `dict`
        `ROOT` : A Path object pointing towards the "root" dir

    Returns:
        new/modified mibig entry as `dict`
    """
    input_message = (
        "================================================\n"
        "Modify the optional information of a MIBiG entry:\n"
        "Enter a number and press enter.\n"
        "Press 'Ctrl+D' to cancel without saving.\n"
        "================================================\n"
        "0) Save and continue\n"
        "1) Gene annotation (partially implemented)\n"
        "2) RiPP annotation (partially implemented)\n"
        "================================================\n"
    )

    while True:
        user_input = input(input_message)

        if user_input == "0":
            break
        elif user_input == "1":
            mibig_entry = get_genes(mibig_entry, ROOT)
            continue
        elif user_input == "2":
            mibig_entry = get_ripp(mibig_entry, ROOT)
            continue
        else:
            print("Invalid input provided")
            continue

    return mibig_entry


def get_changelog(
    mibig_entry: Dict,
    args: argparse.Namespace,
    CURATION_ROUND: str,
    ROOT: Path,
) -> Dict:
    """Collect information for a minimal MIBiG entry.

    Parameters:
        `mibig_entry` : MIBiG entry as `dict`
        `args` : arguments provided by user
        `CURATION_ROUND` : `str` of mibig curation round version to use in changelog
        `ROOT` : `Path` object indicating "root" directory of script

    Returns:
        new/modified mibig entry as `dict`
    """
    curator_id = get_curators(ROOT)[args.curator]

    changelog_object = Changelog(curator_id, CURATION_ROUND)

    if mibig_entry["changelog"][-1]["version"] != CURATION_ROUND:
        return changelog_object.create_new_entry_changelog(mibig_entry)
    else:
        return changelog_object.append_last_entry_changelog(mibig_entry)


def export_mibig_entry(mibig_entry: Dict, path_existing: Path, ROOT: Path) -> None:
    """Concatenate information and create/modify a MIBiG entry.

    Parameters:
        `mibig_dict` : `dict` containing the mibig entry information
        `path_existing` : `Path` object indicating location of optional existing entry
        `ROOT` : `Path` object indicating "root" directory of script

    Returns
        None
    """
    write_mibig = WriteMibig(mibig_entry)

    write_mibig.test_duplicate_entries(ROOT)

    validation_mibig(write_mibig, ROOT)

    if path_existing is None:
        new_entry_path = (
            ROOT.joinpath("mibig_next_ver")
            .joinpath(write_mibig.export_dict["cluster"]["mibig_accession"])
            .with_suffix(".json")
        )
        write_mibig.export_to_json(new_entry_path)
        print("New entry validated and successfully created.")
        write_mibig.append_to_csv_existing(ROOT)
    else:
        write_mibig.export_to_json(path_existing)
        print("Existing entry validated and successfully modified.")

    return


def main() -> None:
    """Entry point of the program.

    For developers:
        This program is designed to ask create/modify MIBiG entries
        by collecting user input, validating it, and writing it in
        the MIBiG JSON format.

        The program is constructed in a modular way: for each "topic"
        of data input, there is a separate class that handles the input
        logic and validates the input. The classes are independent of
        each other so that they can be added and modified as needed.
        There are some exceptions:
            - The Base class, which handles general methods and constants
                which are inherited by the other classes
            - The MibigEntry class, handles the minimum information
                necessary for a MIBiG entry, and which is always executed
            - The Changelog class, which creates/updates the changelog of
                the entry, and which is always executed
            - The WriteMibig class, which converts data into the JSON
                format, and writes the entry to disk
        The script is started in main, where separate functions initialize
        the classes and handle the returned data.
        Further, there are some general functions that are organized
        in the `aux` folder.

        The general data handling concept is as follows: an entry is loaded or
        newly created by the MibigEntry class and handliong of essential data
        takes place. MibigEntry then creates a dict called `mibig_entry`
        which is further modified by the downstream optional classes.
        Once this is finished, the Changelog class adds/modifies the
        changelog, and the WriteMibig class performs final checks and
        writes/exports the file.

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

    mibig_entry = get_optional_data(mibig_entry, ROOT)

    mibig_entry = get_changelog(mibig_entry, args, CURATION_ROUND, ROOT)

    export_mibig_entry(mibig_entry, path_existing, ROOT)


if __name__ == "__main__":
    main()
