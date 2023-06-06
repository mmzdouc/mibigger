#!/usr/bin/env python3

# add class keeping info and reading/writing
# add function initializing class, calling the methods, returning the data
# have a while loop for collecting the data and incrementing

import argparse
from Bio import Entrez
from pathlib import Path
import re
from typing import Dict, Self, List
from urllib.error import HTTPError


from mibig_input_script.aux.error_messages import error_empty_input
from mibig_input_script.aux.error_messages import error_invalid_input
from mibig_input_script.aux.error_messages import error_invalid_char
from mibig_input_script.aux.error_messages import error_var_message


class Minimal:
    """Collects data for MIBiG minimal entry


    TBA
    """

    def __init__(self: Self):
        self.existing_entry = None
        self.biosynth_class: List[str] = []
        self.compound: List[str] = []
        self.accession: str = None  #
        self.completeness = None
        self.start_coord = None  #
        self.end_coord = None  #
        self.publications = None  #
        self.mibig_accession = None
        self.organism_name = None
        self.ncbi_tax_id = None
        self.minimal = None
        self.status = None
        self.evidence: List[str] = []  #

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
                f"{self.mibig_accession}\n"
                "================================================"
            )
            print(message)
        else:
            message = (
                f"================================================\n"
                f"You are CREATING a new MIBiG entry with the ID\n"
                f"{self.mibig_accession}"
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
                # EXPAND
                # function to check if all necessary values are there and to add missing ones
                break
            elif user_input in options:
                options[user_input]()
            else:
                error_invalid_input("option")

    def get_biosynth_class(self: Self) -> None:
        """Get the biosynthetic class of BGC and test if valid

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

        if not len(user_input) == 0:
            biosynth_class = set()
            for selection in user_input:
                if selection in options:
                    biosynth_class.add(options[selection])
                else:
                    error_invalid_input("biosynthetic class selected")
                    return
            self.biosynth_class = list(biosynth_class)
            return
        else:
            error_empty_input("biosynthetic class")
            return

    def get_compound_name(self: Self) -> None:
        """Get the compound name(s)
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

        if not len(user_input) == 0:
            compound = set()
            for i in user_input:
                compound.add(i)
            self.compound = list(compound)
            return
        else:
            error_empty_input("compound name")
            return

    def get_ncbi_data(self: Self) -> None:
        """Get the NCBI accession number, test if valid, get auxiliary info
        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None

        Notes:
            Test suite originally assembled by Barbara Terlouw
            Strictest testing for NCBI Accession number
            Strain ID and TaxID are less problematic - can be added
                by user if they cannot be retrieved automatically
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
            error_empty_input("NCBI Accession number")
            return
        elif any(char in input_accession for char in illegal_chars):
            for char in illegal_chars:
                if char in input_accession:
                    error_invalid_char(char, "NCBI accession number")
                    return
        elif any(input_accession.startswith(chars) for chars in ["GCF_", "GCA_"]):
            error_var_message(f"{input_accession} is an assembly ID")
            return
        elif input_accession.startswith("SRR"):
            error_var_message(
                f"{input_accession} is a SRA record, please enter an assembled contig"
            )
            return
        elif input_accession.startswith("PRJ"):
            error_var_message(
                f"{input_accession} is a bioproject id, please enter a GenBank record"
            )
            return
        elif any(input_accession.startswith(chars) for chars in ["WP_", "YP_"]):
            error_var_message(f"{input_accession} is a protein ID")
            return
        elif input_accession.split(".")[0].endswith("000000"):
            error_var_message(
                f"{input_accession} is a WGS record, please supply the actual contig"
            )
            return
        elif len(input_accession) < 5:
            error_var_message("Accession number too short")
            return
        else:
            try:
                record = Entrez.read(
                    Entrez.efetch(db="nuccore", id=input_accession, retmode="xml")
                )
                self.accession = input_accession
            except HTTPError:
                error_invalid_input("NCBI Accession number")
                return

            try:
                self.organism_name = record[0]["GBSeq_source"]
            except KeyError:
                error_var_message("Organism name not found, is the Accession correct?")

            try:
                tax_record = Entrez.read(
                    Entrez.esearch(db="taxonomy", term=self.organism_name)
                )
                if tax_record["IdList"]:
                    self.ncbi_tax_id = tax_record["IdList"][0]
                else:
                    error_var_message("No NCBI taxonomy ID found")
            except HTTPError:
                error_invalid_input("organism name")

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
            error_invalid_input("start coordinate")
            return

        try:
            input_end = int(input(input_msg_end_coord).replace(" ", ""))
        except ValueError:
            error_invalid_input("end coordinate")
            return

        if input_start >= input_end:
            error_var_message(
                "The end coordinate cannot lie before the start coordinate"
            )
            return
        else:
            self.start_coord = input_start
            self.end_coord = input_end

        return

    def get_organism_data(self: Self) -> None:
        """Get organism and ncbi taxid

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

        if input_strain != "":
            self.organism_name = input_strain
        else:
            error_empty_input("organism name")

        try:
            self.ncbi_tax_id = int(input(input_msg_taxid).replace(" ", ""))
        except ValueError:
            error_invalid_input("NCBI Taxonomy ID")

        return

    def get_evidence(self: Self) -> None:
        """Get evidence for BGC

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

        if not len(user_input) == 0:
            evidence = set()
            for selection in user_input:
                if selection in options:
                    evidence.add(options[selection])
                else:
                    error_invalid_input("evidence provided")
                    return
            self.evidence = list(evidence)
            return
        else:
            error_empty_input("evidence")
            return

    def get_reference(self: Self) -> None:
        """Get evidence for BGC

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None

        Notes

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
            error_empty_input("publication/reference")
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
                    error_invalid_input("DOI found")
                    return
            elif selection == "2":
                input_pmid = input(input_msg_pmid).replace(" ", "")
                if match := re.search(regex_pattern["pmid"], input_pmid):
                    references.append("".join(["pubmed:", match.group(0)]))
                else:
                    error_invalid_input("Pubmed ID found")
                    return
            elif selection == "3":
                input_patent = input(input_msg_patent).replace(" ", "")
                if match := re.search(regex_pattern["patent"], input_patent):
                    references.append("".join(["patent:", match.group(0)]))
                else:
                    error_invalid_input("Patent found")
                    return
            elif selection == "4":
                input_url = input(input_msg_url).replace(" ", "")
                if match := re.search(regex_pattern["url"], input_url):
                    references.append("".join(["url:", match.group(0)]))
            else:
                error_invalid_input("URL found")
                return

        self.publications = references
        return

    def export_attributes_to_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A `dict` to store in MibigEntry() class
        """
        pass


######## CLASS END


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

    return minimal.export_attributes_to_dict()
