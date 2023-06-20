"""
Module for gene annotation.

Module to get information about the genes in the BGC.
"""

from copy import deepcopy
from typing import Dict, Self, List
import re

from mibig_input_script.classes.class_base import BaseClass


class Genes(BaseClass):
    """Module organizing functions to add gene annotations to MIBiG entry.

    Attributes:
        mibig_dict (Dict) : holding the existing mibig entry

    Methods:
        export_dict(self: Self) -> Dict
        initialize_genes_entry(self: Self) -> None
        get_gene_info(self: Self) -> None
        get_annotations(self: Self) -> None
        remove_gene_annotation(self: Self, invar_menu_entries: int) -> None
        get_gene_annot_entry(self: Self, index: int | str) -> None
        write_gene_annot_entry(
            self: Self,
            index: int,
            input_id: str,
            input_name: str,
            input_product: str,
            input_tailoring: List,
            input_comments: str,
        ) -> None
        get_gene_annot_id(self: Self) -> str | bool
        get_gene_annot_name(self: Self) -> str | bool | None
        get_gene_annot_product(self: Self) -> str | None
        get_gene_annot_tailoring(self: Self) -> List | None
        get_gene_annot_comments(self: Self) -> str | None

    Notes:
        Can be expanded for further annotation entries
            - `extra_genes` , `operons` in `get_gene_info()`
            - `domains`, `functions`, `mut_pheno`, `publications`
                in `get_annotations()`
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
        }

        while True:
            try:
                gene_ann = len(self.mibig_dict["cluster"]["genes"]["annotations"])
            except KeyError:
                gene_ann = "None"

            input_message = (
                "================================================\n"
                "Add gene information:\n"
                "Enter a number and press enter.\n"
                "Press 'Ctrl+D' to cancel without saving.\n"
                "================================================\n"
                "0) Save and continue\n"
                f"1) Add/modify gene annotations (currently: '{gene_ann}')\n"
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
        if input_product == False:
            return
        else:
            pass

        input_tailoring = self.get_gene_annot_tailoring()
        if input_tailoring == False:
            return
        else:
            pass

        input_comments = self.get_gene_annot_comments()
        if input_comments == False:
            return
        else:
            pass

        if isinstance(index, str):
            self.mibig_dict["cluster"]["genes"]["annotations"].append({})
            index = -1
            self.write_gene_annot_entry(
                index,
                input_id,
                input_name,
                input_product,
                input_tailoring,
                input_comments,
            )
        elif isinstance(index, int):
            self.write_gene_annot_entry(
                index,
                input_id,
                input_name,
                input_product,
                input_tailoring,
                input_comments,
            )
        else:
            pass

        return

    def write_gene_annot_entry(
        self: Self,
        index: int,
        input_id: str,
        input_name: str,
        input_product: str,
        input_tailoring: List,
        input_comments: str,
    ) -> None:
        """Write gene annotation entry to self.mibig_dict.

        Parameters:
            `self` : The instance of class Gene.

        Returns:
            None
        """
        if input_id is not None:
            self.mibig_dict["cluster"]["genes"]["annotations"][index]["id"] = input_id
        else:
            pass

        if input_name is not None:
            self.mibig_dict["cluster"]["genes"]["annotations"][index][
                "name"
            ] = input_name
        else:
            pass

        if input_product is not None:
            self.mibig_dict["cluster"]["genes"]["annotations"][index][
                "product"
            ] = input_product
        else:
            pass

        if input_tailoring is not None:
            self.mibig_dict["cluster"]["genes"]["annotations"][index][
                "tailoring"
            ] = input_tailoring
        else:
            pass

        if input_comments is not None:
            self.mibig_dict["cluster"]["genes"]["annotations"][index][
                "comments"
            ] = input_comments
        else:
            pass

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
        elif any(re.match(regexp, input_id) for regexp in self.cont_invalid_gene_names):
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
        elif any(
            re.match(regexp, input_name) for regexp in self.cont_invalid_gene_names
        ):
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

    def get_gene_annot_tailoring(self: Self) -> List | None:
        """Get a list of tailoring reactions associated to translated gene.

        Parameters:
            `self` : The instance of class Gene.

        Returns:
            Return a list of tailoring reactions or None (skip)
        """
        input_msg_tailoring = (
            "================================================\n"
            "Enter the tailoring reactions associated to translated gene:\n"
            "To SKIP, press enter.\n"
            "To specify multiple entries, separate them with a TAB character.\n"
            "================================================\n"
        )
        input_tailoring = input(input_msg_tailoring)
        input_tailoring = list(filter(None, input_tailoring.split("\t")))
        if len(input_tailoring) == 0:
            self.message_formatted("Empty input-value: SKIP")
            return None
        else:
            input_tailoring_set = set()
            for entry in input_tailoring:
                input_tailoring_set.add(entry)
            return list(input_tailoring_set)

    def get_gene_annot_comments(self: Self) -> str | None:
        """Get comment regarding gene annotation.

        Parameters:
            `self` : The instance of class Gene.

        Returns:
            Return a str of the comment or None (skip)
        """
        input_msg_comment = (
            "================================================\n"
            "Enter a comment regarding the gene annotation (or SKIP by pressing enter):\n"
            "================================================\n"
        )
        input_comment = input(input_msg_comment)
        if input_comment == "":
            self.message_formatted("Empty input value - SKIP")
            return None
        else:
            return input_comment
