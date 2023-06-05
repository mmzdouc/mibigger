#!/usr/bin/env python3

# add class keeping info and reading/writing
# add function initializing class, calling the methods, returning the data
# have a while loop for collecting the data and incrementing

import argparse
from pathlib import Path
from typing import Dict, Self, List


class Minimal:
    """Collects data for MIBiG minimal entry


    TBA
    """

    def __init__(self: Self):
        self.existing_entry = None
        self.biosynth_class: List[str] = []
        self.compound = None  #
        self.accession = None  #
        self.completeness = None
        self.start_coord = None  #
        self.end_coord = None  #
        self.publications = None  #
        self.mibig_accession = None
        self.organism_name = None
        self.ncbi_tax_id = None
        self.minimal = None
        self.status = None
        self.evidence = None  #

    def load_existing(self: Self, existing: Dict) -> None:
        """Load the data ofa nexisting entry for later manipulation.

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
        self.minimal = existing["cluster"]["minimal"]  # not needed?
        self.status = existing["cluster"]["status"]  # not needed?

        try:
            self.evidence = existing["cluster"]["loci"]["evidence"]
        except KeyError:
            self.evidence = None

    def get_new_mibig_accession(
        self: Self, args: argparse.Namespace, ROOT: Path
    ) -> None:
        """Generate a non-existing temporary MIBiG ID for the new entry.

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
                f"================================================\n"
                f"You are MODIFYING the existing MIBiG entry:\n"
                f"{self.mibig_accession}"
            )
            print(message)
        else:
            message = (
                f"================================================\n"
                f"You are CREATING a new MIBiG entry with the ID\n"
                f"{self.mibig_accession}"
            )
            print(message)

        start_message = (
            "================================================\n"
            "This menu allows to modify the minimum infor-\n"
            "mation of a MIBiG entry. To modify information,\n"
            "type its corresponding number and press enter.\n"
            "Press 'Ctrl+D' to cancel without saving.\n"
            "================================================"
        )
        print(start_message)

        options = {"1": self.get_biosynth_class}

        while True:
            input_message = (
                f"0) Save and continue\n"
                f"1) Biosynthetic class(es) (currently: {self.biosynth_class})\n"
                f"2) Compound name(s) (currently: {self.compound})\n"
                f"3) NCBI accession number (currently: {self.accession})\n"
                f"4) BGC start/end coordinates (currently: {self.start_coord}:{self.end_coord})\n"
                f"5) BGC evidence (currently: {self.evidence})\n"
                f"6) Publication/reference (currently: {self.publications})\n"
                f"================================================\n"
            )
            user_input = input(input_message)

            if user_input == "0":
                # function to check if all necessary values are there and to add missing ones
                break
            elif user_input in options:
                options[user_input]()
            else:
                error_message = (
                    "================================================\n"
                    "Please enter a valid number.\n"
                    "================================================\n"
                )
                print(error_message)

    def get_biosynth_class(self: Self) -> None:
        """Get the biosynthetic class of BGC and test if valid

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        allowed_bgc = [
            "Alkaloid",
            "Polyketide",
            "RiPP",
            "NRP",
            "Saccharide",
            "Terpene",
            "Other",
        ]
        error_message = (
            "================================================\n"
            "Please enter a valid BGC.\n"
            "================================================\n"
        )

        input_message = (
            f"Type in the biosynthetic class of the BGC.\n"
            f"To specify multiple classes, separate entries with whitespace.\n"
            f"Allowed entries are:\n"
            f"------------------------------------------------\n"
            f"{' '.join([i for i in allowed_bgc]) }\n"
            f"------------------------------------------------\n"
        )
        input_raw = input(input_message)
        user_input = input_raw.split()

        if not len(user_input) == 0:
            self.biosynth_class = set()
            for i in user_input:
                if i in allowed_bgc:
                    self.biosynth_class.add(i)
                else:
                    print(error_message)
                    return
            self.biosynth_class = list(self.biosynth_class)
            return
        else:
            print(error_message)
            return

        # Set a flag for "modified" -> write new entry
        # Else, skip writing - abort

    # funciton idea: while loop for both reading and user input the same -> writes to same attributes
    # esearch function for taxid and organism name (via BioPython)

    def create_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A `dict` to store in MibigEntry() class
        """
        pass


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
