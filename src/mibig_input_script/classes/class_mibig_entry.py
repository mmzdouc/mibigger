"""
Module to load/crate mibig entry dict.

For new wntries, creates a new mibig_entry dict.
For existing entries, creates a deep copy of the original entry and
modifies it directly in place. Returns a mibig entry dict that is consequently
modified by downstream modules.
"""

import argparse
from copy import deepcopy
from Bio import Entrez
from pathlib import Path
import re
from typing import Dict, Self, List, Tuple
from urllib.error import HTTPError

from mibig_input_script.classes.class_base import BaseClass


class MibigEntry(BaseClass):
    """Collect data for MIBiG minimal entry.

    Attributes:
        mibig_accession (str | None): existing/new accession number
        mibig_dict (Dict | None): existing/new data for MIBiG entry

    Methods:
        get_new_mibig_accession(
            self: Self, args: argparse.Namespace, ROOT: Path
            ) -> None
        create_new_entry(self: Self, curator: str, CURATION_ROUND: str) -> None
        load_existing_entry(self: Self, existing: Dict) -> None
        get_input(self: Self) -> None
        get_biosynth_class(self: Self) -> None
        get_compound_name(self: Self) -> None
        get_accession_data(self: Self) -> None
        get_organism_data(self: Self) -> None
        get_evidence(self: Self) -> None
        assign_reference(self: Self) -> None
        test_presence_data(self: Self) -> bool
        set_flags(self: Self) -> None
        export_attributes_to_dict(self: Self) -> Dict

    Note:
        Deepcopy required to prevent implicit changing of original dict "existing"
    """

    def __init__(self: Self):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None
        """
        self.mibig_accession: str | None = None
        self.mibig_dict: Dict | None = None

    def get_new_mibig_accession(
        self: Self, args: argparse.Namespace, ROOT: Path
    ) -> None:
        """Generate a non-existing temporary MIBiG ID for the new entry.

        Parameters:
            `self` : The instance of class MibigEntry.
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

        self.mibig_accession = "".join([str(args.curator), str(counter)])
        return

    def create_new_entry(
        self: Self,
        curator: str,
        CURATION_ROUND: str,
    ) -> None:
        """Create a new base MIBiG entry.

        Parameters:
            `self` : The instance of class MibigEntry.
            `curator` : name of curator creating new entry
            `CURATION_ROUND` : `str` of mibig curation round version to use in changelog

        Returns:
            None
        """

        message = (
            f"================================================\n"
            f"You are CREATING a new MIBiG entry with the ID\n"
            f"{self.mibig_accession}\n"
            "================================================"
        )
        print(message)

        self.mibig_dict = {
            "changelog": [
                {
                    "comments": [],
                    "contributors": [curator],
                    "version": CURATION_ROUND,
                }
            ],
            "cluster": {
                "biosyn_class": [],
                "compounds": [],
                "loci": {
                    "accession": "None",
                    "completeness": "None",
                    "start_coord": "None",
                    "end_coord": "None",
                    "evidence": [],
                },
                "mibig_accession": self.mibig_accession,
                "minimal": "None",
                "ncbi_tax_id": "None",
                "organism_name": "None",
                "publications": [],
                "status": "None",
            },
        }
        return

    def load_existing_entry(self: Self, existing: Dict) -> None:
        """Load the data of an existing entry for later manipulation.

        Parameters:
            `self` : The instance of class MibigEntry.
            `existing` : The existing mibig entry

        Returns:
            None
        """
        self.mibig_dict = deepcopy(existing)
        self.mibig_accession = self.mibig_dict["cluster"]["mibig_accession"]

        message = (
            f"================================================\n"
            f"You are MODIFYING the existing MIBiG entry:\n"
            f"{self.mibig_accession}\n"
            "================================================"
        )
        print(message)

        return

    def get_input(self: Self) -> None:
        """Handle methods for user input and data validation.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None
        """
        options = {
            "1": self.get_biosynth_class,
            "2": self.get_compound_name,
            "3": self.get_accession_data,
            "4": self.get_organism_data,
            "5": self.get_evidence,
            "6": self.assign_reference,
        }

        while True:
            # put into a new function "print_selection_menu"

            biosyn_class = self.mibig_dict["cluster"]["biosyn_class"]
            compounds = [
                self.mibig_dict["cluster"]["compounds"][i]["compound"]
                for i in range(len(self.mibig_dict["cluster"]["compounds"]))
            ]
            accession = self.mibig_dict["cluster"]["loci"]["accession"]
            start_coord = self.mibig_dict["cluster"]["loci"]["start_coord"]
            end_coord = self.mibig_dict["cluster"]["loci"]["end_coord"]
            organism_name = self.mibig_dict["cluster"]["organism_name"]
            ncbi_tax_id = self.mibig_dict["cluster"]["ncbi_tax_id"]
            publications = self.mibig_dict["cluster"]["publications"]
            try:
                evidence = self.mibig_dict["cluster"]["loci"]["evidence"]
            except KeyError:
                evidence = "None"

            input_message = (
                "================================================\n"
                "Modify the minimum information of a MIBiG entry:\n"
                "Enter a number and press enter.\n"
                "Press 'Ctrl+D' to cancel without saving.\n"
                "================================================\n"
                "0) Save and continue\n"
                f"1) Biosynthetic class(es) (currently: {biosyn_class})\n"
                f"2) Compound name(s) (currently: {compounds}):\n"
                f"3) NCBI Accession number, start/end coordinates (currently: {accession}, {start_coord}:{end_coord})\n"
                f"4) Organism name, NCBI Taxonomy ID (currently: {organism_name}, {ncbi_tax_id})\n"
                f"5) BGC evidence (currently: {evidence})\n"
                f"6) Publication/reference (currently: {publications})\n"
                f"================================================\n"
            )
            user_input = input(input_message)

            if user_input == "0":
                if self.test_presence_data():
                    self.set_flags()
                    break
                else:
                    self.error_message_formatted(
                        "Not all required information provided"
                    )
                    continue
            elif user_input in options:
                options[user_input]()
                continue
            else:
                self.error_message_formatted("Invalid input provided")
                continue
        return

    def get_biosynth_class(self: Self) -> None:
        """Get the biosynthetic class of BGC and test if valid.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None
        """
        options = {
            "1": "Alkaloid",
            "2": "Polyketide",
            "3": "RiPP",
            "4": "NRP",
            "5": "Saccharide",
            "6": "Terpene",
            "7": "Other",
        }

        input_message = (
            "================================================\n"
            "Separate multiple entries with a TAB character.\n"
            "================================================\n"
            "1) Alkaloid\n"
            "2) Polyketide\n"
            "3) RiPP\n"
            "4) NRP\n"
            "5) Saccharide\n"
            "6) Terpene\n"
            "7) Other\n"
            "================================================\n"
        )

        input_raw = input(input_message)
        user_input = list(filter(None, input_raw.split("\t")))

        if len(user_input) == 0:
            self.error_message_formatted("Empty input value")
            return
        else:
            biosynth_class = set()
            for selection in user_input:
                if selection in options:
                    biosynth_class.add(options[selection])
                else:
                    self.error_message_formatted("Invalid input provided")
                    return
            self.mibig_dict["cluster"]["biosyn_class"] = list(biosynth_class)
            return

    def get_compound_name(self: Self) -> None:
        """Get the compound name(s).

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None
        """
        input_message = (
            "================================================\n"
            "Enter the compound names associated with the BGC.\n"
            "To specify multiple names, separate entries with a TAB character.\n"
            "================================================\n"
        )
        input_raw = input(input_message)
        user_input = list(filter(None, input_raw.split("\t")))

        if len(user_input) == 0:
            self.error_message_formatted("Empty input value")
            return
        else:
            compounds = set()
            for entry in user_input:
                compounds.add(entry)

            compounds = list(compounds)
            self.mibig_dict["cluster"]["compounds"] = []
            for entry in range(len(compounds)):
                self.mibig_dict["cluster"]["compounds"].append(
                    {"compound": compounds[entry]}
                )
            return

    def get_ncbi_accession(self: Self) -> Dict:
        """Get the NCBI accession number, test if valid.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            A `dict` of variables assigned

        Notes:
            Test suite originally assembled by Barbara Terlouw
        """
        input_msg_accession = (
            "================================================\n"
            "Enter the NCBI Accession number (ESearch retrieval takes a few seconds).\n"
            "================================================\n"
        )

        accession_regex = r"^([A-Za-z0-9_]{3,}\.\d)|(MIBIG\.BGC\d{7}\.\d)$"

        illegal_chars = [
            ",",
            "-",
            "|",
            "/",
        ]

        ncbi_data = {
            "ncbi_accession": "None",
            "ncbi_organism_name": "None",
            "ncbi_tax_id": "None",
        }

        input_accession = input(input_msg_accession).replace(" ", "")

        if input_accession == "":
            self.error_message_formatted("Empty input value")
            return ncbi_data
        elif any(char in input_accession for char in illegal_chars):
            for char in illegal_chars:
                if char in input_accession:
                    self.error_message_formatted(
                        f"Illegal character '{char}' in NCBI accession number"
                    )
                    return ncbi_data
        elif any(input_accession.startswith(chars) for chars in ["GCF_", "GCA_"]):
            self.error_message_formatted(f"{input_accession} is an assembly ID")
            return ncbi_data
        elif input_accession.startswith("SRR"):
            self.error_message_formatted(
                f"{input_accession} is a SRA record, please enter an assembled contig"
            )
            return ncbi_data
        elif input_accession.startswith("PRJ"):
            self.error_message_formatted(
                f"{input_accession} is a bioproject id, please enter a GenBank record"
            )
            return ncbi_data
        elif any(input_accession.startswith(chars) for chars in ["WP_", "YP_"]):
            self.error_message_formatted(f"{input_accession} is a protein ID")
            return ncbi_data
        elif input_accession.split(".")[0].endswith("000000"):
            self.error_message_formatted(
                f"{input_accession} is a WGS record, please supply the actual contig"
            )
            return ncbi_data
        elif not re.search(accession_regex, input_accession):
            self.error_message_formatted(
                f"{input_accession} did not match the expected format"
            )
            return ncbi_data
        elif len(input_accession) < 5:
            self.error_message_formatted("Accession number too short")
            return ncbi_data
        else:
            try:
                ncbi_record = Entrez.read(
                    Entrez.efetch(db="nuccore", id=input_accession, retmode="xml")
                )
                ncbi_data["ncbi_accession"] = input_accession
            except HTTPError:
                self.error_message_formatted("Invalid input provided")
                return ncbi_data

            try:
                ncbi_data["ncbi_organism_name"] = ncbi_record[0]["GBSeq_source"]
            except KeyError:
                self.error_message_formatted(
                    "Organism name not found in Accession. Please add manually."
                )
                return ncbi_data

            try:
                tax_record = Entrez.read(
                    Entrez.esearch(db="taxonomy", term=ncbi_data["ncbi_organism_name"])
                )
                if tax_record["IdList"]:
                    ncbi_data["ncbi_tax_id"] = str(tax_record["IdList"][0])
                else:
                    self.error_message_formatted("No NCBI taxonomy ID found")
                    return ncbi_data
            except HTTPError:
                self.error_message_formatted(
                    "Search for NCBI Taxonomy ID via organism name failed"
                )
                return ncbi_data

        return ncbi_data

    def get_coordinates(self: Self) -> Dict:
        """Get the coordinates, test if valid.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            A `dict` of variables assigned
        """
        input_msg_start_coord = (
            "================================================\n"
            "Enter the start coordinate.\n"
            "================================================\n"
        )

        input_msg_end_coord = (
            "================================================\n"
            "Enter the end coordinate.\n"
            "================================================\n"
        )

        coordinates = {"start": "None", "end": "None"}

        try:
            input_start = int(input(input_msg_start_coord).replace(" ", ""))
        except ValueError:
            self.error_message_formatted(
                "Invalid input, start coordinate must be provided"
            )
            return coordinates

        try:
            input_end = int(input(input_msg_end_coord).replace(" ", ""))
        except ValueError:
            self.error_message_formatted(
                "Invalid input, end coordinate must be provided"
            )
            return coordinates

        if input_start >= input_end:
            self.error_message_formatted(
                "The end coordinate cannot lie before the start coordinate"
            )
            return coordinates
        else:
            coordinates["start"] = input_start
            coordinates["end"] = input_end
            return coordinates

    def get_accession_data(self: Self) -> None:
        """Get the NCBI accession number, test if valid, get auxiliary info.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None
        """
        ncbi_data = self.get_ncbi_accession()

        if ncbi_data["ncbi_accession"] != "None":
            self.mibig_dict["cluster"]["loci"]["accession"] = ncbi_data[
                "ncbi_accession"
            ]
            self.mibig_dict["cluster"]["organism_name"] = ncbi_data[
                "ncbi_organism_name"
            ]
            self.mibig_dict["cluster"]["ncbi_tax_id"] = ncbi_data["ncbi_tax_id"]
        else:
            return

        coordinates = self.get_coordinates()

        if coordinates["start"] != "None" and coordinates["end"] != "None":
            self.mibig_dict["cluster"]["loci"]["start_coord"] = coordinates["start"]
            self.mibig_dict["cluster"]["loci"]["end_coord"] = coordinates["end"]
        else:
            return

        return

    def get_organism_data(self: Self) -> None:
        """Get organism and ncbi taxid from user.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None

        Notes:
            Deliberately no testing performed to allow people to correct
            automatically retrieved entries
        """
        input_msg_strain = (
            "================================================\n"
            "Enter the organism name (including strain identifier).\n"
            "================================================\n"
        )

        input_msg_taxid = (
            "================================================\n"
            "Enter the NCBI Taxonomy ID number.\n"
            "================================================\n"
        )

        input_strain = input(input_msg_strain)

        if input_strain == "":
            self.error_message_formatted("Empty input value")
        else:
            self.mibig_dict["cluster"]["organism_name"] = input_strain

        try:
            self.mibig_dict["cluster"]["ncbi_tax_id"] = str(
                int(input(input_msg_taxid).replace(" ", ""))
            )
        except ValueError:
            self.error_message_formatted("Invalid input provided")

        return

    def get_evidence(self: Self) -> None:
        """Get evidence for BGC-compound connection.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        options = {
            "1": "Gene expression correlated with compound production",
            "2": "Knock-out studies",
            "3": "Enzymatic assays",
            "4": "Heterologous expression",
            "5": "In vitro expression",
        }

        input_msg_evidence = (
            "================================================\n"
            "Separate multiple entries with a TAB character.\n"
            "================================================\n"
            "1) Gene expression correlated with compound production\n"
            "2) Knock-out studies\n"
            "3) Enzymatic assays\n"
            "4) Heterologous expression\n"
            "5) In vitro expression\n"
            "================================================\n"
        )

        input_raw = input(input_msg_evidence)
        user_input = list(filter(None, input_raw.split("\t")))

        if len(user_input) == 0:
            self.error_message_formatted("Empty input value")
            return
        else:
            evidence = set()
            for selection in user_input:
                if selection in options:
                    evidence.add(options[selection])
                else:
                    self.error_message_formatted("Invalid input provided")
                    return
            self.mibig_dict["cluster"]["loci"]["evidence"] = list(evidence)
            return

    def assign_reference(self: Self) -> None:
        """Assigns publication reference to BGC.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        reference = self.get_reference()

        if reference is not None:
            self.mibig_dict["cluster"]["publications"] = reference
            return
        else:
            return

    def test_presence_data(self: Self) -> bool:
        """Test if all required attributes are present.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            `bool` to indicate if ready to save entry
        """
        attributes_strs = [
            self.mibig_dict["cluster"]["loci"]["accession"],
            self.mibig_dict["cluster"]["loci"]["start_coord"],
            self.mibig_dict["cluster"]["loci"]["end_coord"],
            self.mibig_dict["cluster"]["mibig_accession"],
            self.mibig_dict["cluster"]["organism_name"],
            self.mibig_dict["cluster"]["ncbi_tax_id"],
        ]

        attributes_lists = [
            self.mibig_dict["cluster"]["biosyn_class"],
            self.mibig_dict["cluster"]["compounds"],
            self.mibig_dict["cluster"]["loci"]["evidence"],
            self.mibig_dict["cluster"]["publications"],
        ]

        if all(variable != "None" for variable in attributes_strs) and (
            all(variable != [] for variable in attributes_lists)
        ):
            return True
        else:
            return False

    def set_flags(self: Self) -> None:
        """Set the appropriate flags to be stored in the MIBiG file.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None

        Notes:
            Set flags for new entries
            For old entries, take presents except for completeness
            (evidence is mandatory in this script)
        """
        if self.mibig_dict["cluster"]["minimal"] == "None":
            self.mibig_dict["cluster"]["minimal"] = True
            self.mibig_dict["cluster"]["status"] = "active"
            self.mibig_dict["cluster"]["loci"]["completeness"] = "complete"
            return
        else:
            self.mibig_dict["cluster"]["loci"]["completeness"] = "complete"
            return

    def export_attributes_to_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A json-compatible dict of MIBiG "cluster" entry.
        """
        return self.mibig_dict
