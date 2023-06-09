"""
Module to handle data to create/modify a MIBiG entry.

Handles dict created by other modules, performs tests to prevent
duplication of entries, writes json files.
"""

from copy import deepcopy
import json
from pathlib import Path
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
        """Placeholder entry to test for duplicate entries when creating a new entry

        Parameters:
            `self` : The instance of class WriteMibig.
            `ROOT` : a `Path` object of the "base" directory

        Returns:
            None
        """
        return
