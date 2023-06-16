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
        get_ncbi_accession(self: Self) -> Dict
        get_coordinates(self: Self) -> Dict
        get_accession_data(self: Self) -> None
        get_organism_data(self: Self) -> None
        get_evidence(self: Self) -> None
        assign_reference(self: Self) -> None
        test_presence_data(self: Self) -> bool
        set_flags(self: Self) -> None
        export_dict(self: Self) -> Dict

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
            `CURATION_ROUND` : `str` of mibig curation round version

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

    def load_existing_entry(self: Self, mibig_entry: Dict) -> None:
        """Load the data of an existing entry for later manipulation.

        Parameters:
            `self` : The instance of class MibigEntry.
            `mibig_entry` : The existing mibig entry

        Returns:
            None
        """
        self.mibig_dict = deepcopy(mibig_entry)
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
                f"2) Compound name(s) (currently: {compounds})\n"
                "3) NCBI Accession number, start:end coordinates "
                f"(currently: '{accession}', '{start_coord}:{end_coord}')\n"
                "4) Organism name, NCBI Taxonomy ID"
                f"(currently: '{organism_name}', '{ncbi_tax_id}')\n"
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
        number_invar_menu_entries = 3

        while True:
            input_message_list = [
                (
                    "================================================\n"
                    "Add new compound or modify existing ones:\n"
                    "Enter a number and press enter.\n"
                    "================================================\n"
                    "0) Save and continue\n"
                    "1) Add a new entry\n"
                    "2) Delete an entry\n"
                )
            ]

            counter = 3

            for entry in range(len(self.mibig_dict["cluster"]["compounds"])):
                input_message_list.append(
                    (
                        f"{counter}) "
                        f"{self.mibig_dict['cluster']['compounds'][entry]['compound']}"
                        f"\n"
                    )
                )
                counter += 1

            input_message_list.append(
                ("================================================\n")
            )

            input_message = "".join([i for i in input_message_list])

            input_raw = input(input_message)

            if input_raw == "":
                self.error_message_formatted("Invalid input provided")
                continue
            elif input_raw == "0":
                break
            elif input_raw == "1":
                self.get_compound_entry(input_raw)
                continue
            elif input_raw == "2":
                self.remove_compound_entry(number_invar_menu_entries)
                continue
            else:
                try:
                    input_raw = int(input_raw) - number_invar_menu_entries
                except ValueError:
                    self.error_message_formatted("Invalid input provided")
                    continue
                if 0 <= input_raw <= len(self.mibig_dict["cluster"]["compounds"]):
                    self.get_compound_entry(input_raw)
                    continue
                else:
                    self.error_message_formatted("Invalid input provided")
                    continue

        return

    def remove_compound_entry(self: Self, number_invar_menu_entries: int) -> None:
        """Remove a compound from the MIBiG entry.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None

        Notes:
        """
        input_msg_number = (
            "================================================\n"
            "Enter the number of compound to remove (a single one):\n"
            "================================================\n"
        )

        input_number = input(input_msg_number)

        try:
            input_number = int(input_number)
        except ValueError:
            self.error_message_formatted("Invalid input value")
            return

        input_number = input_number - number_invar_menu_entries

        if 0 <= input_number <= len(self.mibig_dict["cluster"]["compounds"]):
            if self.ask_proceed_input():
                self.mibig_dict["cluster"]["compounds"].pop(input_number)
                return
            else:
                return
        else:
            self.error_message_formatted("Entry not found")
            return

    def get_compound_entry(self: Self, index: int | str) -> None:
        """Get information on compound.

        Parameters:
            `self` : The instance of class RiPP.
            `str` to create a new entry (append), `int` to replace existing.

        Returns:
            None

        Notes:
            If `index` is `str`, then a new entry is written (appended).
            If `index` is `int`, an existing entry is modified (overwritten).
        """
        input_msg_name = (
            "================================================\n"
            "Enter the compound name (a single one):\n"
            "================================================\n"
        )

        input_msg_synonym = (
            "================================================\n"
            "Enter synonym name(s) (to SKIP, press enter):\n"
            "To specify multiple names, separate entries with a TAB character.\n"
            "================================================\n"
        )

        input_msg_structure = (
            "================================================\n"
            "Enter a single compound structure as SMILES string:\n"
            "(Please create separate entries for congeners)\n"
            "================================================\n"
        )

        input_msg_evidence = (
            "================================================\n"
            "Enter the structural evidence for the chemical structure:\n"
            "(Enter one or more numbers separated by a TAB character and press enter)\n"
            "================================================\n"
            "1) X-ray crystallography\n"
            "2) Nuclear Magnetic Resonance (NMR) spectroscopy\n"
            "3) (Tandem) mass spectrometry (MS/MS)\n"
            "4) Prediction/other\n"
            "================================================\n"
        )

        evidence_options = {
            "1": "X-ray",
            "2": "NMR",
            "3": "MS/MS",
            "4": "Other",
        }

        input_msg_activity = [
            (
                "================================================\n"
                "Enter the biological activities (to SKIP, press enter):\n"
                "(Enter one or more activities separated by a TAB character and press enter)\n"
                "================================================\n"
            )
        ]
        for option in self.const_allowed_bioactiv:
            input_msg_activity.append(f"- {option}\n")
        input_msg_activity.append("================================================\n")
        input_msg_activity = "".join([i for i in input_msg_activity])

        compound_regexp = r"^[a-zA-Zα-ωΑ-Ω0-9\[\]\'()/&,. +-]+$"
        smiles_regexp = r"^[\[\]a-zA-Z0-9@()=\\/\\#+.%*-]+$"

        input_name = input(input_msg_name)
        if input_name == "":
            self.error_message_formatted("Empty input value")
            return
        elif not re.match(compound_regexp, input_name):
            self.error_message_formatted(f"Invalid compound name '{input_name}'")
            return
        else:
            pass

        input_synonym = input(input_msg_synonym)
        input_synonym = list(filter(None, input_synonym.split("\t")))
        if len(input_synonym) == 0:
            self.message_formatted("Empty input value - SKIP")
            input_synonym = None
        else:
            input_synonym_set = set()
            for entry in input_synonym:
                if not re.match(compound_regexp, entry):
                    self.error_message_formatted(f"Invalid compound name '{entry}'")
                    return
                else:
                    input_synonym_set.add(entry)
            input_synonym = list(input_synonym_set)

        input_structure = input(input_msg_structure)
        if input_structure == "":
            self.error_message_formatted("Empty input value")
            return
        elif not re.match(smiles_regexp, input_structure):
            self.error_message_formatted("Invalid SMILES string. Please check again'")
            return
        else:
            pass

        input_evidence = input(input_msg_evidence)
        input_evidence = list(filter(None, input_evidence.split("\t")))
        if len(input_evidence) == 0:
            self.error_message_formatted("Empty input value")
            return
        else:
            input_evidence_set = set()
            for entry in input_evidence:
                try:
                    int(entry)
                except ValueError:
                    self.error_message_formatted(f"Input '{entry}' is not a number")
                    return
                if int(entry) not in [1, 2, 3, 4]:
                    self.error_message_formatted(f"Input '{entry}' not valid")
                    return
                else:
                    input_evidence_set.add(evidence_options[entry])
            input_evidence = list(input_evidence_set)

        input_activity = input(input_msg_activity)
        input_activity = list(filter(None, input_activity.split("\t")))
        if len(input_activity) == 0:
            self.message_formatted("Empty input value - SKIP")
            input_activity = None
        else:
            input_activity_set = set()
            for entry in input_activity:
                if entry not in self.const_allowed_bioactiv:
                    self.error_message_formatted(f"Invalid bioactivity '{entry}'")
                    return
                else:
                    input_activity_set.add(entry)
            input_activity = list(input_activity_set)

        # add target
        if isinstance(index, str):
            self.mibig_dict["cluster"]["compounds"].append(
                {
                    "compound": input_name,
                    "chem_struct": input_structure,
                    "evidence": input_evidence,
                }
            )
            if input_synonym is not None:
                self.mibig_dict["cluster"]["compounds"][-1][
                    "chem_synonyms"
                ] = input_synonym
            if input_activity is not None:
                self.mibig_dict["cluster"]["compounds"][-1]["chem_acts"] = []
                for i in range(len(input_activity)):
                    self.mibig_dict["cluster"]["compounds"][-1]["chem_acts"].append(
                        {"activity": input_activity[i]}
                    )

        elif isinstance(index, int):
            self.mibig_dict["cluster"]["compounds"][index] = {
                "compound": input_name,
                "chem_struct": input_structure,
                "evidence": input_evidence,
            }
            if input_synonym is not None:
                self.mibig_dict["cluster"]["compounds"][index][
                    "chem_synonyms"
                ] = input_synonym
            if input_activity is not None:
                self.mibig_dict["cluster"]["compounds"][index]["chem_acts"] = []
                for i in range(len(input_activity)):
                    self.mibig_dict["cluster"]["compounds"][index]["chem_acts"].append(
                        {"activity": input_activity[i]}
                    )
        else:
            pass

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
            "Enter the NCBI Accession number "
            "(ESearch retrieval takes a few seconds).\n"
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
        elif re.search(r"_", input_accession):
            self.error_message_formatted(
                (
                    f"{input_accession} is a RefSeq record.\n"
                    "Is is recommended to not enter RefSeq records into MIBiG.\n"
                    "Please consider using a GenBank accession number"
                )
            )
            if self.ask_proceed_input():
                pass
            else:
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

    def export_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A json-compatible dict of MIBiG entry
        """
        return deepcopy(self.mibig_dict)
