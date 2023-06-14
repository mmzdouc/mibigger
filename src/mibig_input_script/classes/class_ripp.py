"""
Module to get information about RiPP.

Module to get and store optional RiPP information. Not fully implemented
yet.
"""
from copy import deepcopy
from typing import Dict, Self
import re

from mibig_input_script.classes.class_base import BaseClass


class Ripp(BaseClass):
    """Module organizing functions to add info about RiPP to MIBiG entry

    Attributes:
        mibig_dict (Dict) : holding the existing mibig entry

    Methods:
        initialize_ripp_entry(self: Self) -> None
        set_flag_minimal_false(self: Self) -> None
        get_input(self: Self) -> None
        get_cyclic(self: Self) -> None
        get_subclass(self: Self) -> None
        get_mod_precursor_pep(self: Self) -> None
        get_precursor_peptide_entry(self: Self, index: int | str) -> None
        export_dict(self: Self) -> Dict

    """

    def __init__(self: Self, mibig_entry: Dict):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class RiPP.
            `mibig_entry` : Existing entry to be modified

        Returns:
            None
        """
        self.mibig_dict = deepcopy(mibig_entry)

    def initialize_ripp_entry(self: Self) -> None:
        """Initialize a minimal ripp entry in dict.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        self.mibig_dict["cluster"]["ripp"] = {}
        return

    def set_flag_minimal_false(self: Self) -> None:
        """Set the flag minimal to false.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        self.mibig_dict["cluster"]["minimal"] = False
        return

    def get_input(self: Self) -> None:
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        options = {
            "1": self.get_cyclic,
            "2": self.get_subclass,
            "3": self.get_mod_precursor_pep,
        }

        while True:
            try:
                cyclic = self.mibig_dict["cluster"]["ripp"]["cyclic"]
            except KeyError:
                cyclic = "None"

            try:
                subclass = self.mibig_dict["cluster"]["ripp"]["subclass"]
            except KeyError:
                subclass = "None"

            try:
                mod_precursor_pep = len(
                    self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"]
                )
            except KeyError:
                mod_precursor_pep = "None"

            input_message = (
                "================================================\n"
                "RiPP annotation (currently not fully implemented):\n"
                "Enter a number and press enter.\n"
                "Press 'Ctrl+D' to cancel without saving.\n"
                "================================================\n"
                "0) Save and continue\n"
                f"1) Cyclic (currently: '{cyclic}')\n"
                f"2) RiPP subclass (currently: '{subclass}')\n"
                f"3) Modified precursor peptides (currently: '{mod_precursor_pep}')\n"
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

    def get_cyclic(self: Self) -> None:
        """Get input if cyclic or not.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        input_message = (
            "================================================\n"
            "Is the RiPP cyclic (head-to-tail cyclized)?\n"
            "Enter a number and press enter.\n"
            "================================================\n"
            "1) Cyclic\n"
            "2) Linear\n"
            "================================================\n"
        )

        while True:
            input_raw = input(input_message)
            if input_raw == "1":
                self.mibig_dict["cluster"]["ripp"]["cyclic"] = True
                break
            elif input_raw == "2":
                self.mibig_dict["cluster"]["ripp"]["cyclic"] = False
                break
            else:
                self.error_message_formatted("Invalid input provided")
                continue

        return

    def get_subclass(self: Self) -> None:
        """Get input which subclass of RiPP it is.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None

        Note:
            May be extended to check for allowed subclasses in the future
        """
        input_message = (
            "================================================\n"
            "Enter the subclass of the RiPP.\n"
            "================================================\n"
        )
        while True:
            input_raw = input(input_message)
            if input_raw == "":
                self.error_message_formatted("Invalid input provided")
                continue
            else:
                self.mibig_dict["cluster"]["ripp"]["subclass"] = input_raw
                break
        return

    def get_mod_precursor_pep(self: Self) -> None:
        """Get input on the modified precursor peptide structure.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        try:
            self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"]
        except KeyError:
            self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"] = []

        while True:
            input_message_list = [
                (
                    "================================================\n"
                    "Add new modified precursor peptides or modify existing ones:\n"
                    "Enter a number and press enter.\n"
                    "================================================\n"
                    "0) Save and continue\n"
                    "1) Add new entry\n"
                )
            ]

            counter = 2
            for entry in range(
                len(self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"])
            ):
                input_message_list.append(
                    (
                        f"{counter}) "
                        f"{self.mibig_dict['cluster']['ripp']['mod_precursor_pep'][entry]['name']}"
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
                self.get_precursor_peptide_entry(input_raw)
                continue
            else:
                try:
                    input_raw = int(input_raw) - 2
                except ValueError:
                    self.error_message_formatted("Invalid input provided")
                    continue
                if input_raw <= len(
                    self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"]
                ):
                    self.get_precursor_peptide_entry(input_raw)
                    continue
                else:
                    self.error_message_formatted("Invalid input provided")
                    continue
        return

    def get_precursor_peptide_entry(self: Self, index: int | str) -> None:
        """Get a new modified precursor peptide

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            `str` to create a new entry (append), `int` to replace existing

        Notes:
            If `index` is `str`, then a new entry is written (appended).
            If `index` is `int`, an existing entry is modified (overwritten).
        """
        input_msg_gene_id = (
            "================================================\n"
            "Enter the gene ID:\n"
            "================================================\n"
        )
        input_msg_peptide = (
            "================================================\n"
            "Enter the modified precursor peptide using the shorthand format:\n"
            "(leader)-core-(follower) (e.g. MR-YYH-)\n"
            "================================================\n"
        )
        input_msg_name = (
            "================================================\n"
            "Enter the name of the final product:\n"
            "================================================\n"
        )

        input_gene_id = input(input_msg_gene_id)

        if input_gene_id == "":
            self.error_message_formatted("Invalid input provided")
            return
        elif not re.match(r"^[^, ]*$", input_gene_id):
            self.error_message_formatted("Input value cannot contain space or comma")
            return
        else:
            pass

        input_peptide = input(input_msg_peptide)

        if input_peptide == "":
            self.error_message_formatted("Invalid input provided")
            return
        elif not re.match(r"^[A-Z]*-.+-[A-Z]*$", input_peptide):
            self.error_message_formatted("Modified precursor peptide in wrong format")
            return
        else:
            pass

        input_name = input(input_msg_name)

        if input_peptide == "":
            self.error_message_formatted("Invalid input provided")
            return
        else:
            pass

        if isinstance(index, str):
            self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"].append(
                {"gene_id": input_gene_id, "peptide": input_peptide, "name": input_name}
            )
        elif isinstance(index, int):
            self.mibig_dict["cluster"]["ripp"]["mod_precursor_pep"][index] = {
                "gene_id": input_gene_id,
                "peptide": input_peptide,
                "name": input_name,
            }
        else:
            pass

        return

    def export_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            A json-compatible dict of MIBiG entry
        """
        return deepcopy(self.mibig_dict)
