"""
Module for gene annotation.

Module to get information about the genes in the BGC.
"""

from copy import deepcopy
from typing import Dict, Self, List
import re

from mibig_input_script.classes.class_base import BaseClass


class Genes(BaseClass):
    """Module organizing functions to add gene annotations to MIBiG entry

    Attributes:
        mibig_dict (Dict) : holding the existing mibig entry

    Methods:
        TBA
    """

    def __init__(self: Self, mibig_entry: Dict):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class Genes.
            `mibig_entry` : Existing entry to be modified

        Returns:
            None
        """
        self.mibig_dict = deepcopy(mibig_entry)

    def export_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict.

        Parameters:
            `self` : The instance of class Genes.

        Returns:
            A json-compatible dict of MIBiG entry
        """
        return deepcopy(self.mibig_dict)

    def initialize_genes_entry(self: Self) -> None:
        """Initialize a minimal genes entry.

        Parameters:
            `self` : The instance of class Genes.

        Returns:
            None
        """
        self.mibig_dict["cluster"]["genes"] = {}
        return

    def get_gene_info(self: Self) -> None:
        """Create gene info menu.

        Parameters:
            `self` : The instance of class Genes.

        Returns:
            None
        """
        options = {
            "1": self.get_annotations,
            # ~ "2": self.get_extra_genes,
            # ~ "3": self.get_operons,
        }

        while True:
            try:
                gene_ann = len(self.mibig_dict["cluster"]["genes"]["annotations"])
            except KeyError:
                gene_ann = "None"

            input_message = (
                "================================================\n"
                "Add gene information (currently not fully implemented):\n"
                "Enter a number and press enter.\n"
                "Press 'Ctrl+D' to cancel without saving.\n"
                "================================================\n"
                "0) Save and continue\n"
                f"1) Add/modify gene annotations (currently: '{gene_ann}')\n"
                # ~ "2) Add extra genes\n"
                # ~ "3) Add operon annotations\n"
                "================================================\n"
            )

            user_input = input(input_message)
            if user_input == "0":
                self.set_flag_minimal_false()
                break
            elif user_input in options:
                options[user_input]()
                continue
            else:
                self.error_message_formatted("Invalid input provided")
                continue
        return

    def get_annotations(self: Self) -> None:
        """Create gene annotation menu.

        Parameters:
            `self` : The instance of class Genes.

        Returns:
            None
        """
        try:
            self.mibig_dict["cluster"]["genes"]["annotations"]
        except KeyError:
            self.mibig_dict["cluster"]["genes"]["annotations"] = []

        invar_menu_entries = 3

        while True:
            input_msg_menu = [
                (
                    "================================================\n"
                    "Add new gene annotations or modify existing ones:\n"
                    "Enter a number and press enter.\n"
                    "================================================\n"
                    "0) Save and continue\n"
                    "1) Add new entry\n"
                    "2) Delete an entry\n"
                )
            ]

            counter = 3
            for entry in range(len(self.mibig_dict["cluster"]["genes"]["annotations"])):
                input_msg_menu.append(
                    (
                        f"{counter}) Gene ID: "
                        f"{self.mibig_dict['cluster']['genes']['annotations'][entry]['id']}"
                        f"\n"
                    )
                )
                counter += 1

            input_msg_menu.append(
                ("================================================\n")
            )
            input_msg_menu = "".join([i for i in input_msg_menu])

            input_selection = input(input_msg_menu)

            if input_selection == "":
                self.error_message_formatted("Invalid input provided")
                continue
            elif input_selection == "0":
                break
            elif input_selection == "1":
                self.get_gene_annot_entry(input_selection)
                continue
            elif input_selection == "2":
                self.remove_gene_annotation(invar_menu_entries)
                continue
            else:
                try:
                    input_selection = int(input_selection) - invar_menu_entries
                except ValueError:
                    self.error_message_formatted("Invalid input provided")
                    continue
                if (
                    0
                    <= input_selection
                    < len(self.mibig_dict["cluster"]["genes"]["annotations"])
                ):
                    self.get_gene_annot_entry(input_selection)
                    continue
                else:
                    self.error_message_formatted("Invalid input provided")
                    continue
        return

    def remove_gene_annotation(self: Self, invar_menu_entries: int) -> None:
        """Remove a gene annotation entry from the MIBiG entry.

        Parameters:
            `self` : The instance of class Genes.

        Returns:
            None
        """
        input_msg_number = (
            "================================================\n"
            "Enter the number of the gene annotation to remove (a single one):\n"
            "================================================\n"
        )
        input_number = input(input_msg_number)

        try:
            input_number = int(input_number)
        except ValueError:
            self.error_message_formatted("Invalid input value")
            return

        input_number = input_number - invar_menu_entries

        if 0 <= input_number < len(self.mibig_dict["cluster"]["genes"]["annotations"]):
            if self.ask_proceed_input():
                self.mibig_dict["cluster"]["genes"]["annotations"].pop(input_number)
                return
            else:
                return
        else:
            self.error_message_formatted("Entry not found")
            return

    def get_gene_annot_entry(self: Self, index: int | str) -> None:
        """Add a gene annotation entry.

        Parameters:
            `self` : The instance of class Genes.
            `str` to create a new entry (append), `int` to replace existing.

        Returns:
            None

        Notes:
            If `index` is `str`, then a new entry is written (appended).
            If `index` is `int`, an existing entry is modified (overwritten).
        """

        input_id = self.get_gene_annot_id()
        if input_id == False:
            return
        else:
            pass

        input_name = self.get_gene_annot_name()
        if input_name == False:
            return
        else:
            pass

        input_product = self.get_gene_annot_product()

        print(input_product)

        # id : str
        # name : str
        # product : str
        # tailoring : List
        # comments : str

        return

    def get_gene_annot_id(self: Self) -> str | bool:
        """Get the protein ID of the translated gene.

        Parameters:
            `self` : The instance of class Genes.

        Returns:
            Return a str of the gene id or
            bool (error)
        """
        input_msg_id = (
            "================================================\n"
            "Enter the protein ID of the translated gene:\n"
            "================================================\n"
        )
        input_id = input(input_msg_id)
        if input_id == "":
            self.error_message_formatted("Protein ID cannot be empty")
            return False
        elif not re.match(r"^[^, ]*$", input_id):
            self.error_message_formatted(f"'{input_id}' is not in the correct format")
            return False
        else:
            return input_id

    def get_gene_annot_name(self: Self) -> str | bool | None:
        """Get the gene name.

        Parameters:
            `self` : The instance of class Gene.

        Returns:
            Return a str of the gene name or
            bool (error) or None (skip)
        """
        input_msg_name = (
            "================================================\n"
            "Enter the gene name (or SKIP by pressing enter):\n"
            "================================================\n"
        )
        input_name = input(input_msg_name)
        if input_name == "":
            self.message_formatted("Empty input value - SKIP")
            return None
        elif not re.match(r"^[^, ]*$", input_name):
            self.error_message_formatted(f"'{input_name}' is not in the correct format")
            return False
        else:
            return input_name

    def get_gene_annot_product(self: Self) -> str | None:
        """Get the gene product description.

        Parameters:
            `self` : The instance of class Gene.

        Returns:
            Return a str of the gene product description or None (skip)
        """
        input_msg_product = (
            "================================================\n"
            "Enter a description of the gene product (or SKIP by pressing enter):\n"
            "================================================\n"
        )
        input_product = input(input_msg_product)
        if input_product == "":
            self.message_formatted("Empty input value - SKIP")
            return None
        else:
            return input_product

    # ~ continue here
