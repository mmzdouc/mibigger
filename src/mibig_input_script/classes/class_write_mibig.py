"""
Module to handle data to create/modify a MIBiG entry

Handles dict created by other modules, performs tests to prevent
duplication of entries, writes json files
"""

from copy import deepcopy
from pathlib import Path
from typing import Dict, Self


class WriteMibig:
    """Class to prepare, test, and write a new/existing MIBiG entry to disk.

    Attributes:
        TBA

    Methods:
        TBA

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
        self.export_dict: Dict = None

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

    def test_duplicate_entries(self: Self, ROOT) -> None:
        """Check existing entries for similar-looking entry"""
        pass

    def return_json_string(self: Self, ROOT) -> None:
        """Validate modified entry before writing"""
        # return a json string to test with validation script of Simon/Kai
        pass

    def export_to_json(self: Self, path_to_write: Path) -> None:
        """Save to json"""
        pass
