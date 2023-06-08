"""
Module to create/modify a MIBiG changelog.

Manages new/existing entries, keeps track of changes and who performed them
For new entries, write a new changelog
For existing entries not modified in current curation round, create new
changelog entry
For existing entries already modified in current curation round, append
to the last changelog entry.
"""

from datetime import datetime
from typing import Set, Dict, List, Self


class Changelog:
    """Collect data for MIBiG minimal entry.

    Attributes:
        existing_entry (Dict | None): optional existing entry to append to.
        curator (str) : Curator modifying entry
        set_curators (Set) : Set of curators modifying entry in curation round
        version (str) : current version
        comments (List) : Comment explaining modification

    Methods:
        error_message_formatted(self: Self, string: str) -> None:
            Print a formatted error message.
        get_comment(self: Self) -> None:
            Get a changelog comment from curator.
        create_new_changelog(self: Self) -> Dict:
            Create a new changelog from scratch.
        create_new_entry_changelog(self: Self, existing: Dict) -> Dict:
            Add a new changelog entry to existing changelog.
        append_last_entry_changelog(self: Self, existing: Dict) -> Dict:
            Append info to last entry of existing changelog.
    """

    def __init__(self: Self, curator: str, CURATION_ROUND: str):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class Minimal.
            `curator` : A `str` indicating the curator working on entry
            `CURATION_ROUND` : A `str` indicating the current curation round

        Returns:
            None
        """
        self.existing_entry: Dict | None = None
        self.curator: str = curator
        self.set_curators: Set = set()
        self.version: str = CURATION_ROUND
        self.comments: List = list()

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

    def get_comment(self: Self) -> None:
        """Get a changelog comment from curator.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        input_msg_comment = (
            "================================================\n"
            "Please enter a brief comment summarizing the changes.\n"
            "================================================\n"
        )

        while True:
            input_comment = str(input(input_msg_comment))
            if input_comment == "":
                self.error_message_formatted("Comment cannot be empty")
                continue
            else:
                now = datetime.now()
                formatted_dt = now.strftime("%d-%m-%Y %H:%M")
                self.comments.append(
                    input_comment + f" ({self.curator} {formatted_dt})"
                )
                break

    def create_new_changelog(self: Self) -> Dict:
        """Create a new changelog from scratch.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            A json-compatible dict of a "changelog" MIBiG entry.
        """
        self.get_comment()

        return {
            "changelog": [
                {
                    "comments": self.comments,
                    "contributors": [self.curator],
                    "version": self.version,
                }
            ]
        }

    def create_new_entry_changelog(self: Self, existing: Dict) -> Dict:
        """Add a new changelog entry to existing changelog.

        Parameters:
            `self` : The instance of class Minimal.
            `existing` : The loaded MIBiG json entry as `dict`.

        Returns:
            A json-compatible dict of a "changelog" MIBiG entry.
        """
        self.existing_entry = existing
        self.get_comment()

        self.existing_entry["changelog"].append(
            {
                "comments": self.comments,
                "contributors": [self.curator],
                "version": self.version,
            }
        )

        return {"changelog": self.existing_entry["changelog"]}

    def append_last_entry_changelog(self: Self, existing: Dict) -> Dict:
        """Append info to last entry of existing changelog.

        Parameters:
            `self` : The instance of class Minimal.
            `existing` : The loaded MIBiG json entry as `dict`.

        Returns:
            A json-compatible dict of a "changelog" MIBiG entry.
        """
        self.existing_entry = existing

        self.comments = existing["changelog"][-1]["comments"]
        self.get_comment()

        self.set_curators.update(existing["changelog"][-1]["contributors"])
        self.set_curators.add(self.curator)

        self.existing_entry["changelog"][-1]["comments"] = self.comments
        self.existing_entry["changelog"][-1]["contributors"] = list(self.set_curators)

        return {"changelog": self.existing_entry["changelog"]}
