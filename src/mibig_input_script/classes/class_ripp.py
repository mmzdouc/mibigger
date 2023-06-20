"""
Module to get information about RiPP.

Module to get and store optional RiPP information.
"""
from copy import deepcopy
from typing import Dict, Self, List
import re

from mibig_input_script.classes.class_base import BaseClass


class Ripp(BaseClass):
    """Module organizing functions to add info about RiPP to MIBiG entry.

    Attributes:
        mibig_dict (Dict) : holding the existing mibig entry

    Methods:
        export_dict(self: Self) -> Dict
        initialize_ripp_entry(self: Self) -> None
        get_input(self: Self) -> None
        get_cyclic(self: Self) -> None
        get_subclass(self: Self) -> None
        get_peptidases(self: Self) -> None
        get_precursor_genes(self: Self) -> None
        remove_precursor_gene_entry(self: Self, invar_menu_entries: int) -> None
        get_precursor_gene_entry(self: Self, index: int | str) -> None
        write_precursor_gene_entry(
            self: Self,
            index : int,
            input_leader : str,
            input_core : List,
            input_follower : str,
            input_gene_id : str,
            input_cleavage : List,
            input_recognition : str,
        ) -> None
        get_precursor_gene_leader(self: Self) -> str | bool | None
        get_precursor_gene_core(self: Self) -> List | bool
        get_precursor_gene_follower(self: Self) -> str | bool | None
        get_precursor_gene_id(self: Self) -> str | bool
        get_precursor_gene_cleavage(self: Self) -> List | bool | None
        get_precursor_recognition(self: Self) -> str | bool | None

    Notes:
        Can be expanded for further annotation entries
            - `crosslinks` in `get_precursor_genes()`
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

    def export_dict(self: Self) -> Dict:
        """Summarize values in json-compatible dict.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            A json-compatible dict of MIBiG entry
        """
        return deepcopy(self.mibig_dict)

    def initialize_ripp_entry(self: Self) -> None:
        """Initialize a minimal ripp entry in dict.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        self.mibig_dict["cluster"]["ripp"] = {}
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
            "3": self.get_peptidases,
            "4": self.get_precursor_genes,
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
                peptidases = self.mibig_dict["cluster"]["ripp"]["peptidases"]
            except KeyError:
                peptidases = "None"

            try:
                precursor_genes = len(
                    self.mibig_dict["cluster"]["ripp"]["precursor_genes"]
                )
            except KeyError:
                precursor_genes = "None"

            input_message = (
                "================================================\n"
                "RiPP annotation (currently not fully implemented):\n"
                "Enter a number and press enter.\n"
                "Press 'Ctrl+D' to cancel without saving.\n"
                "================================================\n"
                "0) Save and continue\n"
                f"1) Cyclic (currently: '{cyclic}')\n"
                f"2) RiPP subclass (currently: '{subclass}')\n"
                f"3) RiPP-cleaving peptidases (currently: '{peptidases}')\n"
                f"4) RiPP precursor genes (currently: '{precursor_genes}')\n"
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
        input_msg_cyclic = (
            "================================================\n"
            "Is the RiPP cyclic (head-to-tail cyclized)?\n"
            "Enter a number and press enter.\n"
            "================================================\n"
            "1) Cyclic\n"
            "2) Linear\n"
            "================================================\n"
        )
        input_cyclic = input(input_msg_cyclic)
        if input_cyclic == "1":
            self.mibig_dict["cluster"]["ripp"]["cyclic"] = True
            return
        elif input_cyclic == "2":
            self.mibig_dict["cluster"]["ripp"]["cyclic"] = False
            return
        else:
            self.error_message_formatted("Invalid input provided")
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
        input_msg_subclass = (
            "================================================\n"
            "Enter the subclass of the RiPP (or SKIP by pressing enter).\n"
            "================================================\n"
        )

        input_subclass = input(input_msg_subclass)
        if input_subclass == "":
            self.message_formatted("Invalid input provided - SKIP")
            return
        else:
            self.mibig_dict["cluster"]["ripp"]["subclass"] = input_subclass
            return

    def get_peptidases(self: Self) -> None:
        """Get peptidase(s) involved in precursor cleavage.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        input_msg_peptidases = (
            "================================================\n"
            "Enter the peptidase(s) involved in precursor cleavage.\n"
            "To SKIP, press enter.\n"
            "To specify multiple entries, separate them with a TAB character.\n"
            "================================================\n"
        )
        input_peptidases = input(input_msg_peptidases)
        input_peptidases = list(filter(None, input_peptidases.split("\t")))

        if len(input_peptidases) == 0:
            self.message_formatted("Empty input value - SKIP")
            return
        else:
            peptidases_set = set()
            for entry in input_peptidases:
                peptidases_set.add(entry)
            self.mibig_dict["cluster"]["ripp"]["peptidases"] = list(peptidases_set)
            return

    def get_precursor_genes(self: Self) -> None:
        """Get input on precursor genes.

        Parameters:
            `self` : The instance of class RiPP.

        Returns:
            None
        """
        try:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"]
        except KeyError:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"] = []

        invar_menu_entries = 3

        while True:
            input_msg_menu = [
                (
                    "================================================\n"
                    "Add new precursor genes or modify existing ones:\n"
                    "Enter a number and press enter.\n"
                    "================================================\n"
                    "0) Save and continue\n"
                    "1) Add new entry\n"
                    "2) Delete an entry\n"
                )
            ]

            counter = 3
            for entry in range(
                len(self.mibig_dict["cluster"]["ripp"]["precursor_genes"])
            ):
                input_msg_menu.append(
                    (
                        f"{counter}) Gene ID: "
                        f"{self.mibig_dict['cluster']['ripp']['precursor_genes'][entry]['gene_id']}"
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
                self.get_precursor_gene_entry(input_selection)
                continue
            elif input_selection == "2":
                self.remove_precursor_gene_entry(invar_menu_entries)
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
                    < len(self.mibig_dict["cluster"]["ripp"]["precursor_genes"])
                ):
                    self.get_precursor_gene_entry(input_selection)
                    continue
                else:
                    self.error_message_formatted("Invalid input provided")
                    continue

        return

    def remove_precursor_gene_entry(self: Self, invar_menu_entries: int) -> None:
        """Remove a precursor gene entry from the MIBiG entry.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            None
        """
        input_msg_number = (
            "================================================\n"
            "Enter the number of precursor gene to remove (a single one):\n"
            "================================================\n"
        )
        input_number = input(input_msg_number)

        try:
            input_number = int(input_number)
        except ValueError:
            self.error_message_formatted("Invalid input value")
            return

        input_number = input_number - invar_menu_entries

        if (
            0
            <= input_number
            < len(self.mibig_dict["cluster"]["ripp"]["precursor_genes"])
        ):
            if self.ask_proceed_input():
                self.mibig_dict["cluster"]["ripp"]["precursor_genes"].pop(input_number)
                return
            else:
                return
        else:
            self.error_message_formatted("Entry not found")
            return

    def get_precursor_gene_entry(self: Self, index: int | str) -> None:
        """Add a precursor gene entry.

        Parameters:
            `self` : The instance of class MibigEntry.
            `str` to create a new entry (append), `int` to replace existing.

        Returns:
            None

        Notes:
            If `index` is `str`, then a new entry is written (appended).
            If `index` is `int`, an existing entry is modified (overwritten).
        """
        input_leader = self.get_precursor_gene_leader()
        if input_leader == False:
            return
        else:
            pass

        input_core = self.get_precursor_gene_core()
        if input_core == False:
            return
        elif input_core is None:
            self.error_message_formatted("Core peptide cannot be empty")
            return
        else:
            pass

        input_follower = self.get_precursor_gene_follower()
        if input_follower == False:
            return
        else:
            pass

        input_gene_id = self.get_precursor_gene_id()
        if input_gene_id == False:
            return
        elif input_gene_id is None:
            self.error_message_formatted("Gene ID cannot be empty")
            return
        else:
            pass

        input_cleavage = self.get_precursor_gene_cleavage()
        if input_cleavage == False:
            return
        else:
            pass

        input_recognition = self.get_precursor_recognition()
        if input_recognition == False:
            return
        else:
            pass

        if isinstance(index, str):
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"].append({})
            index = -1
            self.write_precursor_gene_entry(
                index,
                input_leader,
                input_core,
                input_follower,
                input_gene_id,
                input_cleavage,
                input_recognition,
            )
        elif isinstance(index, int):
            self.write_precursor_gene_entry(
                index,
                input_leader,
                input_core,
                input_follower,
                input_gene_id,
                input_cleavage,
                input_recognition,
            )
        else:
            pass

        return

    def write_precursor_gene_entry(
        self: Self,
        index: int,
        input_leader: str,
        input_core: List,
        input_follower: str,
        input_gene_id: str,
        input_cleavage: List,
        input_recognition: str,
    ) -> None:
        """Write precursor gene entry data to self.mibig_dict.

        Parameters:
            `self` : The instance of class MibigEntry.
            `index` : The index int of precursor gene to write data to
            `input_leader` : A str of the leader peptide sequence
            `input_core` : A List of core sequences
            `input_follower` : A str of the follower peptide sequence
            `input_gene_id` : A str of the gene ID of the precursor gene
            `input_cleavage` : A List of cleavage recognition sites
            `input_recognition` : A str of the enzyme recognition motif

        Returns:
            None
        """
        if input_leader is not None:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"][index][
                "leader_sequence"
            ] = input_leader
        else:
            pass

        if input_core is not None:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"][index][
                "core_sequence"
            ] = input_core
        else:
            pass

        if input_follower is not None:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"][index][
                "follower_sequence"
            ] = input_follower
        else:
            pass

        if input_gene_id is not None:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"][index][
                "gene_id"
            ] = input_gene_id
        else:
            pass

        if input_cleavage is not None:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"][index][
                "cleavage_recogn_site"
            ] = input_cleavage
        else:
            pass

        if input_recognition is not None:
            self.mibig_dict["cluster"]["ripp"]["precursor_genes"][index][
                "recognition_motif"
            ] = input_recognition
        else:
            pass

        return

    def get_precursor_gene_leader(self: Self) -> str | bool | None:
        """Get the precursor gene leader sequence.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            Return a str of the precursor gene leader sequence or
            bool (error) or None
        """
        input_msg_leader = (
            "================================================\n"
            "Enter the leader peptide sequence (or SKIP by pressing enter):\n"
            "Only enter proteinogenic amino acids (upper-case letters).\n"
            "================================================\n"
        )
        input_leader = input(input_msg_leader)
        if input_leader == "":
            self.message_formatted("Empty input value - SKIP")
            return None
        elif not re.match(self.const_proteinogenic_aa_regexp, input_leader):
            self.error_message_formatted(
                f"'{input_leader}' is not in single-letter proteinogenic amino acid code"
            )
            return False
        else:
            return input_leader

    def get_precursor_gene_core(self: Self) -> List | bool:
        """Get a list of precursor gene core sequences.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            Return the list of precursor gene core sequences or bool (error)
        """
        input_msg_core = (
            "================================================\n"
            "Enter the core peptide sequence(s):\n"
            "Only enter proteinogenic amino acids (upper-case letters).\n"
            "To specify multiple ones, separate entries with a TAB character.\n"
            "================================================\n"
        )
        input_core = input(input_msg_core)
        input_core = list(filter(None, input_core.split("\t")))
        if len(input_core) == 0:
            self.error_message_formatted("Core peptide sequences cannot be empty")
            return False
        else:
            input_core_set = set()
            for entry in input_core:
                if not re.match(self.const_proteinogenic_aa_regexp, entry):
                    self.error_message_formatted(
                        f"'{entry}' is not in single-letter proteinogenic amino acid code"
                    )
                    return False
                else:
                    input_core_set.add(entry)
            return list(input_core_set)

    def get_precursor_gene_follower(self: Self) -> str | bool | None:
        """Get the precursor gene follower sequence.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            Return a str of the precursor gene follower sequence or
            bool (error) or None (skip)
        """
        input_msg_follower = (
            "================================================\n"
            "Enter the follower peptide sequence (or SKIP by pressing enter):\n"
            "Only enter proteinogenic amino acids (upper-case letters).\n"
            "================================================\n"
        )
        input_follower = input(input_msg_follower)
        if input_follower == "":
            self.message_formatted("Empty input value - SKIP")
            return None
        elif not re.match(self.const_proteinogenic_aa_regexp, input_follower):
            self.error_message_formatted(
                f"'{input_follower}' is not in single-letter proteinogenic amino acid code"
            )
            return False
        else:
            return input_follower

    def get_precursor_gene_id(self: Self) -> str | bool:
        """Get the precursor gene follower sequence.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            Return a str of the precursor gene id or
            bool (error)
        """
        input_msg_gene_id = (
            "================================================\n"
            "Enter the precursor gene ID:\n"
            "Only enter proteinogenic amino acids (upper-case letters).\n"
            "================================================\n"
        )
        input_precursor_gene_id = input(input_msg_gene_id)
        if input_precursor_gene_id == "":
            self.error_message_formatted("The precursor gene ID cannot be empty!")
            return False
        elif not re.match(r"^[^, ]*$", input_precursor_gene_id):
            self.error_message_formatted(
                f"'{input_precursor_gene_id}' has invalid format"
            )
            return False
        else:
            return input_precursor_gene_id

    def get_precursor_gene_cleavage(self: Self) -> List | bool | None:
        """Get a list of cleavage recognition sites.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            Return the list of cleavage recognition sites or bool (error) or None (skip)
        """
        input_msg_cleavage = (
            "================================================\n"
            "Enter cleavage recognition site(s) (to SKIP, press enter):\n"
            "Only enter proteinogenic amino acids (upper-case letters).\n"
            "To specify multiple ones, separate entries with a TAB character.\n"
            "================================================\n"
        )
        input_cleavage = input(input_msg_cleavage)
        input_cleavage = list(filter(None, input_cleavage.split("\t")))
        if len(input_cleavage) == 0:
            self.message_formatted("Empty input value - SKIP")
            return None
        else:
            input_cleavage_set = set()
            for entry in input_cleavage:
                if not re.match(self.const_proteinogenic_aa_regexp, entry):
                    self.error_message_formatted(
                        f"'{entry}' is not in single-letter proteinogenic amino acid code"
                    )
                    return False
                else:
                    input_cleavage_set.add(entry)
            return list(input_cleavage_set)

    def get_precursor_recognition(self: Self) -> str | bool | None:
        """Get the precursor gene recognition sequence.

        Parameters:
            `self` : The instance of class MibigEntry.

        Returns:
            Return a str of the precursor gene recognition sequence or
            bool (error) or None (skip)
        """
        input_msg_recognition = (
            "================================================\n"
            "Enter the recognition motif in the leader peptide (or SKIP by pressing enter):\n"
            "Only enter proteinogenic amino acids (upper-case letters).\n"
            "================================================\n"
        )
        input_recognition = input(input_msg_recognition)
        if input_recognition == "":
            self.message_formatted("Empty input value - SKIP")
            return None
        elif not re.match(self.const_proteinogenic_aa_regexp, input_recognition):
            self.error_message_formatted(
                f"'{input_recognition}' is not in single-letter proteinogenic amino acid code"
            )
            return False
        else:
            return input_recognition
