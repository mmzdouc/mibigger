"""Base class to organize commonly used functions.

Several functions in the individual classes are (almost) identical,
e.g. error printing, testing for correct input formats. To keep code
DRY, these functions are organized in the base class, which the other
classes inherit.
"""

import json
from pathlib import Path
from typing import Self, List, Dict
import re


class BaseClass:
    """Organizes commonly used functions for inheritance by other classes.

    To adhere to DRY principles, commonly used functions e.g. for
    input testing are organized in this class for inheritance by
    more specialized downstream classes.

    All class attributes start with `const_` to distinguish them (defined
    in single place to be able to adjust centrally)

    Class attributes:
        const_allowed_bioactiv (List) : allowed bioactivities
        const_compound_regexp (str) : regexp allowed compound name format
        const_smiles_regexp (str) : regexp SMILES format
        const_compound_evidence (str) : allowed compound evidence categories

    Methods:
        message_formatted(self: Self, string: str) -> None
            Print a formatted message
        error_message_formatted(self: Self, string: str) -> None
            Print a formatted error message.
        ask_proceed_input(self: Self) -> None
            Ask user to proceed/stop.
        get_reference(self: Self) -> List | None
            Get reference and test for correct format.
        read_known_values_from_json(self: Self, location: Path) -> Dict
            Read a json file containing a list of allowed entries.
        write_known_values_to_json(self: Self, location: Path, values: Dict) -> None
            Write the dict containing the known values to file.


    """

    cont_invalid_gene_names = [
        "no_accession",
        "unknown",
    ]

    const_allowed_bioactiv: List = [
        "adhesion",
        "anthelmintic",
        "antialgal",
        "antibacterial",
        "antibacterial (Gram-negative)",
        "antibacterial (Gram-positive)",
        "anticancer",
        "anticoccidial",
        "antifungal",
        "antiinflammatory",
        "antimalarial",
        "antineoplastic",
        "antioomycete",
        "antioxidant",
        "antiparasidal",
        "antiplasmodial",
        "antiproliferative",
        "antiprotozoal",
        "antitubulin",
        "antitumor",
        "antiviral",
        "biofilm",
        "carbon storage",
        "cell differentiation",
        "cell envelope",
        "cell protectant",
        "cell wall",
        "cold stress",
        "cyst formation",
        "cytotoxic",
        "denitrificative",
        "dermatotoxic",
        "DNA-interfering",
        "emulsifier",
        "enterotoxic",
        "exopolysaccharide",
        "extracelluar capsule",
        "flavor",
        "fluorescent",
        "hemolytic",
        "hepatotoxic",
        "herbicidal",
        "immunomodulatory",
        "immunosuppressive",
        "inducer",
        "inhibitor",
        "insecticidal",
        "ionophore",
        "iron reducing",
        "irritant",
        "morphogen",
        "neuroprotective",
        "neurotoxic",
        "nitrogen reduction",
        "odorous metabolite",
        "osmolytic",
        "phytotoxic",
        "pigment",
        "predation",
        "radical scavenging",
        "regulatory",
        "siderophore",
        "signalling",
        "sodium channel blocking",
        "surfactant",
        "swarming motility",
        "toxic",
        "tumor promoter",
        "UV protective",
        "vesicant",
        "virulence factor",
    ]

    const_compound_regexp = r"^[a-zA-Zα-ωΑ-Ω0-9\[\]\'()/&,. +-]+$"
    const_smiles_regexp = r"^[\[\]a-zA-Z0-9@()=\\/\\#+.%*-]+$"
    const_doi_regexp = r"10\.\d{4,9}/[-\._;()/:a-zA-Z0-9]+"
    const_pubmed_id_regexp = r"\d+"
    const_url_regexp = r"https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*)"
    const_ncbi_accession_regexp = r"^([A-Za-z0-9_]{3,}\.\d)|(MIBIG\.BGC\d{7}\.\d)$"
    const_proteinogenic_aa_regexp = r"^[A-Z]+$"

    const_ncbi_acc_illegal_chars = [
        ",",
        "-",
        "|",
        "/",
    ]

    const_compound_evidence = {
        "1": "X-ray",
        "2": "NMR",
        "3": "MS/MS",
        "4": "Other",
    }

    const_biosyn_classes = {
        "1": "Alkaloid",
        "2": "Polyketide",
        "3": "RiPP",
        "4": "NRP",
        "5": "Saccharide",
        "6": "Terpene",
        "7": "Other",
    }

    const_bgc_evidence = {
        "1": "Correlation of genomic and metabolomic data",
        "2": "Homology-based prediction",
        "3": "Gene expression correlated with compound production",
        "4": "Knock-out studies",
        "5": "Enzymatic assays",
        "6": "Heterologous expression",
        "7": "In vitro expression",
    }

    def __init__(self: Self):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        pass

    def message_formatted(self: Self, string: str) -> None:
        """Print a formatted message.

        Parameters:
            `self` : The instance of class Minimal.
            `string` : input to customize message

        Returns:
            None
        """
        message = (
            "================================================\n"
            f"{string}.\n"
            "================================================\n"
        )
        print(message)
        return

    def error_message_formatted(self: Self, string: str) -> None:
        """Print a formatted error message.

        Parameters:
            `self` : The instance of class Minimal.
            `string` : input to customize message

        Returns:
            None
        """
        error_message = (
            "++++++++++++++++++++++++++++++++++++++++++++++++\n"
            f"ERROR: {string}.\n"
            "++++++++++++++++++++++++++++++++++++++++++++++++"
        )
        print(error_message)
        return

    def ask_proceed_input(self: Self) -> bool:
        """Ask user to proceed/stop.

        Parameters:
            `self` : The instance of class BaseClass.

        Returns:
            A bool from user input
        """
        alert_message = (
            "================================================\n"
            "Do you want to proceed (yes/no)?\n"
            "================================================\n"
        )
        while True:
            user_input = input(alert_message)

            if user_input == "no":
                return False
            elif user_input == "yes":
                return True
            else:
                print("Please type 'yes' or 'no'.")
                continue

    def set_flag_minimal_false(self: Self) -> None:
        """Set the flag minimal to false.

        Parameters:
            `self` : The instance of class Base.

        Returns:
            None
        """
        self.mibig_dict["cluster"]["minimal"] = False
        return

    def read_known_values_from_json(self: Self, location: Path) -> Dict:
        """Read a json file containing a list of allowed entries.

        Parameters:
            `self` : The instance of class RiPP.
            `location` : The Path of the file to be loaded.

        Returns:
            A dict containing the list of allowed entries
        """
        with open(location) as f:
            return json.load(f)

    def write_known_values_to_json(self: Self, location: Path, values: Dict) -> None:
        """Write the dict containing the known values to file.

        Parameters:
            `self` : The instance of class RiPP.
            `location` : The Path of the file to be stored.
            `values` : A Dict containing the new values.

        Returns:
            None
        """
        with open(location, "w", encoding="utf-8") as f:
            json.dump(values, f, indent=4, ensure_ascii=False)
            return
