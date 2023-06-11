"""
Module to handle data to create/modify a MIBiG entry.

Handles dict created by other modules, performs tests to prevent
duplication of entries, writes json files.
"""

from copy import deepcopy
import json
from pathlib import Path
import pandas as pd
from typing import Dict, Self


class WriteMibig:
    """Class to prepare, test, and write a new/existing MIBiG entry to disk.

    Attributes:
        export_dict (Dict) : stores the dictionary for export

    Methods:
        error_message_formatted(self: Self, string: str) -> None:
            Print a formatted error message
        concatenate_dicts(self: Self, mibig_dict: Dict, changelog: Dict) -> None:
            Concatenate individual dicts for export
    Note:
        Deepcopy required to prevent implicit changing of original dict "existing"
    """

    def __init__(self: Self):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class WriteMibig.

        Returns:
            None
        """
        self.export_dict: Dict | None = None

    def error_message_formatted(self: Self, string: str) -> None:
        """Print a formatted error message.

        Parameters:
            `self` : The instance of class Changelog.
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

    def alert_message_formatted(self: Self, string: str, var: str) -> None:
        """Print a formatted alert message.

        Parameters:
            `self` : The instance of class Changelog.
            string : input to customize message


        Returns:
            None
        """
        alert_message = (
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
            f"ALERT: {string}.\n"
            f"{var}\n"
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        )
        print(alert_message)
        return

    def ask_proceed_input(self: Self) -> None:
        """Ask user to proceed/stop.

        Parameters:
            `self` : The instance of class Changelog.

        Returns:
            None
        """
        alert_message = (
            "================================================\n"
            "Do you want to proceed (yes/no)?\n"
            "================================================\n"
        )
        while True:
            user_input = input(alert_message)

            if user_input == "no":
                quit()
            elif user_input == "yes":
                break
            else:
                print("Please type in yes or no.")

        return

    def concatenate_dicts(self: Self, mibig_dict: Dict, changelog: Dict) -> None:
        """Concatenate individual dicts for export.

        Parameters:
            `self` : The instance of class WriteMibig.
            `mibig_dict` : a `dict` containing the minimal MIBiG info
            `changelog` : a `dict` containing the MIBiG changelog

        Returns:
            None

        Notes:
            Once script goes beyond minimal input (e.g. ripp_dict):
            If ripp_dict not None:
                concatenate to export_dict
        """
        self.export_dict = {
            "changelog": deepcopy(changelog["changelog"]),
            "cluster": deepcopy(mibig_dict["cluster"]),
        }

    def return_json_string(self: Self) -> str:
        """Returns a json string.

        Parameters:
            `self` : The instance of class WriteMibig.

        Returns:
            Returns a json string
        """
        return json.dumps(self.export_dict, indent=4, ensure_ascii=False)

    def export_to_json(self: Self, path_to_write: Path) -> None:
        """Save MIBiG entry to json file using the specified path.

        Parameters:
            `self` : The instance of class WriteMibig.
            `path_to_write` : a `Path` object pointing towards file location

        Returns:
            None
        """
        with open(path_to_write, "w", encoding="utf-8") as outfile:
            json.dump(self.export_dict, outfile, indent=4, ensure_ascii=False)

    def test_duplicate_entries(self: Self, ROOT) -> None:
        """Test for duplicate entries when creating a new entry.

        Parameters:
            `self` : The instance of class WriteMibig.
            `ROOT` : a `Path` object of the "base" directory

        Returns:
            None
        """
        existing = "existing_mibig_entries.csv"
        df = pd.read_csv(ROOT.joinpath(existing))

        current_mibig_acc = self.export_dict["cluster"]["mibig_accession"]
        ncbi_acc = self.export_dict["cluster"]["loci"]["accession"]
        compounds = [
            self.export_dict["cluster"]["compounds"][i]["compound"]
            for i in range(len(self.export_dict["cluster"]["compounds"]))
        ]

        if not (
            matches := df.loc[df["accession"].str.contains(ncbi_acc) == True]
        ).empty:
            if current_mibig_acc not in matches["mibig_accession"].to_list():
                self.alert_message_formatted(
                    f"Similar existing accession number in '{existing}' found", matches
                )
                self.ask_proceed_input()
            else:
                return
        else:
            pass

        for compound in compounds:
            if not (
                matches := df.loc[df["compounds"].str.contains(compound) == True]
            ).empty:
                if current_mibig_acc not in matches["mibig_accession"].to_list():
                    self.alert_message_formatted(
                        f"Similar existing compound name in '{existing}' found", matches
                    )
                    self.ask_proceed_input()
                else:
                    return
            else:
                return

    def append_to_csv_existing(self: Self, ROOT) -> None:
        """Append info of new entry to existing_mibig_entries.csv

        Parameters:
            `self` : The instance of class WriteMibig.
            `ROOT` : a `Path` object of the "base" directory

        Returns:
            None
        """
        existing = "existing_mibig_entries.csv"
        df = pd.read_csv(ROOT.joinpath(existing))

        current_mibig_acc = self.export_dict["cluster"]["mibig_accession"]
        ncbi_acc = self.export_dict["cluster"]["loci"]["accession"]
        compounds = "|".join(
            [
                self.export_dict["cluster"]["compounds"][i]["compound"]
                for i in range(len(self.export_dict["cluster"]["compounds"]))
            ]
        )

        new_row = {
            "mibig_accession": current_mibig_acc,
            "compounds": compounds,
            "accession": ncbi_acc,
        }

        df2 = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

        df2.to_csv(existing, index=False)

        return
