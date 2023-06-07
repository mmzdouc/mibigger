"""
Module to create/modify a MIBiG changelog

Manages new/existing entries, keeps track of changes and who performed them
"""

# ~ For new entries, create a new changelog with the tag "next".
# ~ For existing entries, if the entry was not modified in this round of
# ~ curation, add a new changelog-entry with the tag "next".
# ~ If a tag with "next" already exists, check if the curator is different
# ~ from the current one, and only if different, write to changelog


# Each update needs a comment + contributor + date
# Contributors should be a set to keep contributors unique
# Whenever an entry is modified (i.e. saved, add another comment)


# have two dict export function: one for new entry, one for existing
# one a real dict export, the other modifies the existing entry in situ (can be called "append")

from typing import Set, Dict, List, Self


class Changelog:
    """Collect data for MIBiG minimal entry.

    Attributes:
        existing_entry (Dict | None): optional existing entry to append to.
        curator (Set | None) : Curator modifying entry in this curation round
        version (str | None) : current version
        comments (str | None) : Comment explaining modification

    Method:
        load_existing
        get_comment
        append_existing
        create_first_entry
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

    def load_existing(self: Self, existing: Dict) -> None:
        """Load the data of an existing entry for later manipulation.

        Parameters:
            `self` : The instance of class Minimal.
            `existing` : The loaded MIBiG json entry as `dict`.

        Returns:
            None
        """
