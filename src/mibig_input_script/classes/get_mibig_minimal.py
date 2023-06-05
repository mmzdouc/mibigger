#!/usr/bin/env python3

# add class keeping info and reading/writing
# add function initializing class, calling the methods, returning the data
# have a while loop for collecting the data and incrementing

import argparse
from pathlib import Path
from typing import Dict, Self


class Minimal:
    """Collects data for MIBiG minimal entry


    TBA
    """

    def __init__(self: Self):
        self.existing_entry = None
        self.biosynth_class = None
        self.compound = None
        self.accession = None
        self.completeness = None
        self.start_coord = None
        self.end_coord = None
        self.publications = None
        self.mibig_accession = None
        self.organism_name = None
        self.ncbi_tax_id = None
        self.minimal = None
        self.status = None

        # need to by try/except
        self.chem_struct = None
        self.evidence = None

    def load_existing(self: Self, existing: Dict) -> None:
        """Load and store an existing MIBiG entry for manipulation

        Parameters:
            `self` : The instance of class Minimal.
            `existing` : The loaded MIBiG json entry as `dict`.

        Returns:
            None
        """
        self.existing_entry = existing

        self.biosynth_class = existing["cluster"]["biosyn_class"]
        self.compound = existing["cluster"]["compounds"][0]["compound"]
        self.accession = existing["cluster"]["loci"]["accession"]
        self.completeness = existing["cluster"]["loci"]["completeness"]
        self.start_coord = existing["cluster"]["loci"]["start_coord"]
        self.end_coord = existing["cluster"]["loci"]["end_coord"]
        self.publications = existing["cluster"]["publications"]
        self.mibig_accession = existing["cluster"]["mibig_accession"]
        self.organism_name = existing["cluster"]["organism_name"]
        self.ncbi_tax_id = existing["cluster"]["ncbi_tax_id"]
        self.minimal = existing["cluster"]["minimal"]
        self.status = existing["cluster"]["status"]

        try:
            self.chem_struct = existing["cluster"]["compounds"][0]["chem_struct"]
        except KeyError:
            self.chem_struct = None

        try:
            self.evidence = existing["cluster"]["compounds"][0]["evidence"]
        except KeyError:
            self.evidence = None

    def get_new_mibig_accession(
        self: Self, args: argparse.Namespace, ROOT: Path
    ) -> None:
        """Create a MIBiG accession number for the new entry.

        Parameters:
            `self` : The instance of class Minimal.
            `args` : arguments provided by user
            `ROOT` : `Path` object indicating "root" directory of script

        Returns:
            None
        """
        counter = 1

        while Path.exists(
            ROOT.joinpath("mibig_next_ver")
            .joinpath("_".join([str(args.curator), str(counter)]))
            .with_suffix(".json")
        ):
            counter += 1

        self.mibig_accession = "_".join([str(args.curator), str(counter)])

    def get_input(self: Self) -> None:
        """While loop to accept user input and display current attributes

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        if self.existing_entry is not None:
            message = (
                f"===========================================\n"
                f"You are MODIFYING the existing MIBiG entry:\n"
                f"{self.mibig_accession}\n"
                f"Press 'Ctrl+D' to cancel without saving.\n"
                f"===========================================\n"
            )
            print(message)
        else:
            message = (
                f"==============================================\n"
                f"You are CREATING a new MIBiG entry with the ID\n"
                f"{self.mibig_accession}\n"
                f"Press 'Ctrl+D' to cancel without saving.\n"
                f"==============================================\n"
            )
            print(message)

        # ~ while True:

        # print options menu (which entry to add/change)
        # this then calls a separate function with input field and input testing and variable assignment
        # then goes back to top until

        # START WHILE LOOP HERE!

        # Set a flag for "modified" -> write new entry
        # Else, skip writing - abort

    def test_placeholder(self: Self, input_value: str) -> bool:
        """Perform test on input value

        Parameters:
            `self` : The instance of class Minimal.
            `input_value` : The input value to test.

        Returns:
            A `bool` to indicate outcome of test
        """
        pass

    def create_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A `dict` to store in MibigEntry() class
        """
        pass

    # funciton idea: while loop for both reading and user input the same -> writes to same attributes
    # esearch function for taxid and organism name (via BioPython)


# redo this function - already read


def get_mibig_minimal(
    existing_mibig: Dict | None, args: argparse.Namespace, ROOT: Path
) -> Dict:
    """Collect information for a minimal MIBiG entry

    Parameters:
        `existing_mibig` : existing MIBiG entry as `dict` or None
        `args` : arguments provided by user
        `ROOT` : `Path` object indicating "root" directory of script

    Returns:
        TBA

    """

    minimal = Minimal()

    if existing_mibig is not None:
        minimal.load_existing(existing_mibig)
        minimal.get_input()
        # while loop for manipulating existing -> the same as below?

    else:
        minimal.get_new_mibig_accession(args, ROOT)
        minimal.get_input()

    # instead of reading all values of json, keep it as it is
    # modify the values in situ if they are needed

    # go over all self-entries that are not None and add them to
    # original entry before writing it out

    # the "minimal" entry is then dynamic and does not need to be
    # specified completely - can be tested by json schema anyway

    # in while loop, use only overwriting

    return minimal.create_dict()
