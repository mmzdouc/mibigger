"""Base class to organize commonly used functions.

Several functions in the individual classes are (almost) identical,
e.g. error printing, testing for correct input formats. To keep code
DRY, these functions are organized in the base class, which the other
classes inherit.
"""

import re
from typing import Self, List


class BaseClass:
    """Organizes commonly used functions for inheritance by other classes.

    To adhere to DRY principles, commonly used functions e.g. for
    input testing are organized in this class for inheritance by
    more specialized downstream classes.

    Attributes:
        None

    Methods

    """

    def __init__(self: Self):
        """Initialize class attributes.

        Parameters:
            `self` : The instance of class Minimal.

        Returns:
            None
        """
        pass

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
            "++++++++++++++++++++++++++++++++++++++++++++++++\n"
        )
        print(error_message)
        return

    def ask_proceed_input(self: Self) -> None:
        """Ask user to proceed/stop.

        Parameters:
            `self` : The instance of class BaseClass.

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
                print("Abort process.")
                quit()
            elif user_input == "yes":
                break
            else:
                print("Please type 'yes' or 'no'.")

        return

    def get_reference(self: Self) -> List | None:
        """Get reference and test for correct format.

        Parameters:
            `self` : The instance of class BaseClass.

        Returns:
            A list of reference entries OR None.
        """
        input_msg_reference = (
            "================================================\n"
            "Choose which reference to add.\n"
            "Separate multiple entries with a TAB character.\n"
            "================================================\n"
            "1) Digital Object Identifier (DOI - strongly preferred).\n"
            "2) Pubmed ID.\n"
            "3) Patent reference.\n"
            "4) URL.\n"
            "================================================\n"
        )
        input_msg_doi = (
            "================================================\n"
            "Enter a DOI.\n"
            "================================================\n"
        )
        input_msg_pmid = (
            "================================================\n"
            "Enter a Pubmed ID.\n"
            "================================================\n"
        )
        input_msg_patent = (
            "================================================\n"
            "Enter a patent reference.\n"
            "================================================\n"
        )
        input_msg_url = (
            "================================================\n"
            "Enter an URL.\n"
            "================================================\n"
        )
        regex_pattern = {
            "doi": r"10\.\d{4,9}/[-\._;()/:a-zA-Z0-9]+",
            "pmid": r"\d+",
            "patent": r".+",
            "url": r"https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*)",
        }

        input_raw = input(input_msg_reference)
        user_input = list(filter(None, input_raw.split("\t")))

        if len(user_input) == 0:
            self.error_message_formatted("Empty input value")
            return None
        else:
            pass

        reference = list()

        for selection in user_input:
            if selection == "1":
                input_doi = input(input_msg_doi).replace(" ", "")
                if match := re.search(regex_pattern["doi"], input_doi):
                    reference.append("".join(["doi:", match.group(0)]))
                else:
                    self.error_message_formatted("DOI has the wrong fromat")
                    return None
            elif selection == "2":
                input_pmid = input(input_msg_pmid).replace(" ", "")
                if match := re.search(regex_pattern["pmid"], input_pmid):
                    reference.append("".join(["pubmed:", match.group(0)]))
                else:
                    self.error_message_formatted("Pubmed ID has the wrong format")
                    return None
            elif selection == "3":
                input_patent = input(input_msg_patent).replace(" ", "")
                if match := re.search(regex_pattern["patent"], input_patent):
                    reference.append("".join(["patent:", match.group(0)]))
                else:
                    self.error_message_formatted("Patent has the wrong format")
                    return None
            elif selection == "4":
                input_url = input(input_msg_url).replace(" ", "")
                if match := re.search(regex_pattern["url"], input_url):
                    reference.append("".join(["url:", match.group(0)]))
                else:
                    self.error_message_formatted("URL has the wrong format")
                    return None
            else:
                self.error_message_formatted("Invalid input provided")
                return None

        return reference
