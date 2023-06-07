"""
Module to read and handle minimum data necessary for a MIBiG entry.

This module provides classes and functions for the reading and handling
of the minimal data necessary for creating a MIBiG entry. Further, it
allows for the manipulation of existing entries. Each type of user input
is organized in a separate function, including input sanitation and checks.
"""


import argparse
from Bio import Entrez
from pathlib import Path
import re
from typing import Dict, Self, List
from urllib.error import HTTPError


class Minimal:
    """Collect data for MIBiG minimal entry.

    Attributes:
        existing_entry (Dict | None): dict of optional existing entry.
        biosynth_class (List | None): List of biosynthetic classes.
        compound (List | None): List of compound name(s).
        accession (str | None): NCBI Accession number.
        start_coord (int | None): Start coordinate for BGC in Accession.
        end_coord (int | None): End coordinate for BGC in Accession.
        publications (List | None): References associated with the BGC.
        mibig_accession (str | None): MIBiG Accession ID of entry.
        organism_name (str | None): Organism strain name.
        ncbi_tax_id (str | None): NCBI taxonomy ID of organism.
        evidence (List | None): Evidence for BGC-compound connection.
        minimal (bool | None): Flag for a minimal entry.
        status (str | None): Flag for status of entry (active/retired).
        completeness (str | None): Flag for completeness of locus information.

    Methods:
        load_existing(self: Self, existing: Dict) -> None:
            Load the data of an existing entry for later manipulation.
        get_new_mibig_accession(
                self: Self, args: argparse.Namespace, ROOT: Path
                ) -> None:
            Generate a non-existing temporary MIBiG ID for the new entry.
        get_input(self: Self) -> None:
            Handle methods for user input and data validation
        get_biosynth_class(self: Self) -> None:
            Get the biosynthetic class of BGC and test if valid
        get_compound_name(self: Self) -> None:
            Get the compound name(s)
        get_ncbi_data(self: Self) -> None:
            Get the NCBI accession number, test if valid, get auxiliary info
        get_organism_data(self: Self) -> None:
            Get organism and ncbi taxid from user
        get_evidence(self: Self) -> None:
            Get evidence for BGC-compound connection
        get_reference(self: Self) -> None:
            Get publication/reference for BGC
        test_presence_attributes(self: Self) -> bool:
            Test if all required attributes are present
        set_flags(self: Self) -> None:
            Set the appropriate flags to be stored in the MIBiG file
        export_attributes_to_dict(self: Self) -> Dict:
            Summarize values in json-compatible dict
    """

    def __init__(self: Self):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        self.existing_entry: Dict | None = None
        self.biosynth_class: List | None = None
        self.compound: List | None = None
        self.accession: str | None = None
        self.start_coord: int | None = None
        self.end_coord: int | None = None
        self.publications: List | None = None
        self.mibig_accession: str | None = None
        self.organism_name: str | None = None
        self.ncbi_tax_id: str | None = None
        self.evidence: List | None = None
        self.minimal: bool | None = None
        self.status: str | None = None
        self.completeness: str | None = None

    def error_message_formatted(self: Self, string: str) -> None:
        """Print a formatted error message.

        Parameters:
            `self` : The instance of class Minimal.
            string : input to customize message

        Returns:
            None
        """
        error_message = (
            "++++++++++++++++++++++++++++++++++++++++++++++++\n"
            f"ERROR: {string}.\n"
            "++++++++++++++++++++++++++++++++++++++++++++++++\n"
        )
        print(error_message)
        return

    def load_existing(self: Self, existing: Dict) -> None:
        """Load the data of an existing entry for later manipulation.

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
        self.start_coord = existing["cluster"]["loci"]["start_coord"]
        self.end_coord = existing["cluster"]["loci"]["end_coord"]
        self.publications = existing["cluster"]["publications"]
        self.mibig_accession = existing["cluster"]["mibig_accession"]
        self.organism_name = existing["cluster"]["organism_name"]
        self.ncbi_tax_id = existing["cluster"]["ncbi_tax_id"]

        try:
            self.evidence = existing["cluster"]["loci"]["evidence"]
        except KeyError:
            self.evidence = None

        self.minimal = existing["cluster"]["minimal"]
        self.status = existing["cluster"]["status"]
        self.completeness = existing["cluster"]["loci"]["completeness"]

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
        """Handle methods for user input and data validation.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        if self.existing_entry is not None:
            message = (
                f"================================================\n"
                f"You are MODIFYING the existing MIBiG entry:\n"
                f"{self.mibig_accession}\n"
                "================================================"
            )
            print(message)
        else:
            message = (
                f"================================================\n"
                f"You are CREATING a new MIBiG entry with the ID\n"
                f"{self.mibig_accession}\n"
                "================================================"
            )
            print(message)

        options = {
            "1": self.get_biosynth_class,
            "2": self.get_compound_name,
            "3": self.get_ncbi_data,
            "4": self.get_organism_data,
            "5": self.get_evidence,
            "6": self.get_reference,
        }

        while True:
            input_message = (
                "================================================\n"
                "Modify the minimum information of a MIBiG entry:\n"
                "Enter a number and press enter.\n"
                "Press 'Ctrl+D' to cancel without saving.\n"
                "================================================\n"
                f"0) Save and continue\n"
                f"1) Biosynthetic class(es) (currently: {self.biosynth_class})\n"
                f"2) Compound name(s) (currently: {self.compound})\n"
                f"3) NCBI Accession number, start/end coordinates\n"
                f"   (currently: {self.accession}, {self.start_coord}:{self.end_coord})\n"
                f"4) Organism name, NCBI Taxonomy ID (currently: {self.organism_name}, {self.ncbi_tax_id})\n"
                f"5) BGC evidence (currently: {self.evidence})\n"
                f"6) Publication/reference (currently: {self.publications})\n"
                f"================================================\n"
            )
            user_input = input(input_message)

            if user_input == "0":
                if self.test_presence_attributes():
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

    def get_biosynth_class(self: Self) -> None:
        """Get the biosynthetic class of BGC and test if valid.

        Parameters:
            `self` : The instance of class Minimal.

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
            self.biosynth_class = list(biosynth_class)
            return

    def get_compound_name(self: Self) -> None:
        """Get the compound name(s).

        Parameters:
            `self` : The instance of class Minimal.

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
            for i in user_input:
                compounds.add(i)
            self.compound = list(compounds)
            return

    def get_ncbi_data(self: Self) -> None:
        """Get the NCBI accession number, test if valid, get auxiliary info.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None

        Notes:
            Test suite originally assembled by Barbara Terlouw
        """
        input_msg_accession = (
            "================================================\n"
            "Enter the NCBI Accession number (ESearch retrieval takes a few seconds).\n"
            "================================================\n"
        )

        illegal_chars = [
            ",",
            "-",
            "|",
            "/",
        ]

        input_accession = input(input_msg_accession).replace(" ", "")

        if input_accession == "":
            self.error_message_formatted("Empty input value")
            return
        elif any(char in input_accession for char in illegal_chars):
            for char in illegal_chars:
                if char in input_accession:
                    self.error_message_formatted(
                        f"Illegal character '{char}' in NCBI accession number"
                    )
                    return
        elif any(input_accession.startswith(chars) for chars in ["GCF_", "GCA_"]):
            self.error_message_formatted(f"{input_accession} is an assembly ID")
            return
        elif input_accession.startswith("SRR"):
            self.error_message_formatted(
                f"{input_accession} is a SRA record, please enter an assembled contig"
            )
            return
        elif input_accession.startswith("PRJ"):
            self.error_message_formatted(
                f"{input_accession} is a bioproject id, please enter a GenBank record"
            )
            return
        elif any(input_accession.startswith(chars) for chars in ["WP_", "YP_"]):
            self.error_message_formatted(f"{input_accession} is a protein ID")
            return
        elif input_accession.split(".")[0].endswith("000000"):
            self.error_message_formatted(
                f"{input_accession} is a WGS record, please supply the actual contig"
            )
            return
        elif len(input_accession) < 5:
            self.error_message_formatted("Accession number too short")
            return
        else:
            try:
                record = Entrez.read(
                    Entrez.efetch(db="nuccore", id=input_accession, retmode="xml")
                )
                self.accession = input_accession
            except HTTPError:
                self.error_message_formatted("Invalid input provided")
                return

            try:
                self.organism_name = record[0]["GBSeq_source"]
            except KeyError:
                self.error_message_formatted(
                    "Organism name not found, is the Accession correct?"
                )

            try:
                tax_record = Entrez.read(
                    Entrez.esearch(db="taxonomy", term=self.organism_name)
                )
                if tax_record["IdList"]:
                    self.ncbi_tax_id = str(tax_record["IdList"][0])
                else:
                    self.error_message_formatted("No NCBI taxonomy ID found")
            except HTTPError:
                self.error_message_formatted("Invalid input provided")

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

        try:
            input_start = int(input(input_msg_start_coord).replace(" ", ""))
        except ValueError:
            self.error_message_formatted("Invalid input provided")
            return

        try:
            input_end = int(input(input_msg_end_coord).replace(" ", ""))
        except ValueError:
            self.error_message_formatted("Invalid input provided")
            return

        if input_start >= input_end:
            self.error_message_formatted(
                "The end coordinate cannot lie before the start coordinate"
            )
            return
        else:
            self.start_coord = input_start
            self.end_coord = input_end

        return

    def get_organism_data(self: Self) -> None:
        """Get organism and ncbi taxid from user.

        Parameters:
            `self` : The instance of class Minimal.

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
            self.organism_name = input_strain

        try:
            self.ncbi_tax_id = str(int(input(input_msg_taxid).replace(" ", "")))
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
            self.evidence = list(evidence)
            return

    def get_reference(self: Self) -> None:
        """Get publication/reference for BGC.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        regex_pattern = {
            "doi": r"10\.\d{4,9}/[-\._;()/:a-zA-Z0-9]+",
            "pmid": r"\d+",
            "patent": r".+",
            "url": r"https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*)",
        }

        input_msg_reference = (
            "================================================\n"
            "Choose which reference to add.\n"
            "Separate multiple entries with a TAB character.\n"
            "================================================\n"
            "1) Digital Object Identifier (DOI - strongly preferred).\n"
            "2) Pubmed ID.\n"
            "3) Patent reference.\n"
            "4) URL.\n"
            "================================================\n"
        )
        input_msg_doi = (
            "================================================\n"
            "Enter a DOI.\n"
            "================================================\n"
        )
        input_msg_pmid = (
            "================================================\n"
            "Enter a Pubmed ID.\n"
            "================================================\n"
        )
        input_msg_patent = (
            "================================================\n"
            "Enter a patent reference.\n"
            "================================================\n"
        )
        input_msg_url = (
            "================================================\n"
            "Enter an URL.\n"
            "================================================\n"
        )

        input_raw = input(input_msg_reference)
        user_input = list(filter(None, input_raw.split("\t")))

        if len(user_input) == 0:
            self.error_message_formatted("Empty input value")
            return
        else:
            pass

        references = list()

        for selection in user_input:
            if selection == "1":
                input_doi = input(input_msg_doi).replace(" ", "")
                if match := re.search(regex_pattern["doi"], input_doi):
                    references.append("".join(["doi:", match.group(0)]))
                else:
                    self.error_message_formatted("DOI has the wrong fromat")
                    return
            elif selection == "2":
                input_pmid = input(input_msg_pmid).replace(" ", "")
                if match := re.search(regex_pattern["pmid"], input_pmid):
                    references.append("".join(["pubmed:", match.group(0)]))
                else:
                    self.error_message_formatted("Pubmed ID has the wrong format")
                    return
            elif selection == "3":
                input_patent = input(input_msg_patent).replace(" ", "")
                if match := re.search(regex_pattern["patent"], input_patent):
                    references.append("".join(["patent:", match.group(0)]))
                else:
                    self.error_message_formatted("Patent has the wrong format")
                    return
            elif selection == "4":
                input_url = input(input_msg_url).replace(" ", "")
                if match := re.search(regex_pattern["url"], input_url):
                    references.append("".join(["url:", match.group(0)]))
                else:
                    self.error_message_formatted("URL has the wrong format")
                    return
            else:
                self.error_message_formatted("Invalid input provided")
                return

        self.publications = references
        return

    def test_presence_attributes(self: Self) -> bool:
        """Test if all required attributes are present.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            `bool` to indicate if ready to save entry
        """
        attributes = [
            self.biosynth_class,
            self.compound,
            self.accession,
            self.start_coord,
            self.end_coord,
            self.publications,
            self.mibig_accession,
            self.organism_name,
            self.ncbi_tax_id,
            self.evidence,
        ]

        if all(variable is not None for variable in attributes):
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
            For a new entry, all three flag are None. "Minimal" and "Status"
            can be set to True and "active", respectively.
            "Completeness" can also be set to "complete" since the locus
            info needs to be complete before the entry can be saved
            (checked in self.test_completeness_minimal())
        """
        if self.minimal is None:
            self.minimal = True
            self.status = "active"
            self.completeness = "complete"
            return
        else:
            return

    def export_attributes_to_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A `dict` to store in MibigEntry() class
        """
        return {
            "cluster": {
                "biosyn_class": self.biosynth_class,
                "compounds": [{"compound": self.compound}],
                "loci": {
                    "accession": self.accession,
                    "completeness": self.completeness,
                    "start_coord": self.start_coord,
                    "end_coord": self.end_coord,
                    "evidence": self.evidence,
                },
                "mibig_accession": self.mibig_accession,
                "minimal": self.minimal,
                "ncbi_tax_id": self.ncbi_tax_id,
                "organism_name": self.organism_name,
                "publications": self.publications,
                "status": self.status,
            }
        }
