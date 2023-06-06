#!/usr/bin/env python3#


def error_empty_input(string: str) -> None:
    """Print an error message due to empty input

    Parameters:
        string : input to customize message

    Returns:
        None
    """
    error_message = (
        "++++++++++++++++++++++++++++++++++++++++++++++++\n"
        f"ERROR: The {string} cannot be empty!\n"
        "++++++++++++++++++++++++++++++++++++++++++++++++\n"
    )

    print(error_message)
    return None


def error_invalid_input(string: str) -> None:
    """Print an error message due to invalid input

    Parameters:
        string : input to customize message

    Returns:
        None
    """
    error_message = (
        "++++++++++++++++++++++++++++++++++++++++++++++++\n"
        f"ERROR: Not a valid {string}!\n"
        "++++++++++++++++++++++++++++++++++++++++++++++++\n"
    )

    print(error_message)
    return None


def error_invalid_char(char: str, attribute: str) -> None:
    """Print an error message due to invalid input

    Parameters:
        char : input to customize message
        attribute : input to customize message

    Returns:
        None
    """
    error_message = (
        "++++++++++++++++++++++++++++++++++++++++++++++++\n"
        f"ERROR: {char} is not a valid character in a {attribute}!\n"
        "++++++++++++++++++++++++++++++++++++++++++++++++\n"
    )

    print(error_message)
    return None


def error_var_message(
    string: str,
) -> None:
    """Print an error message due to invalid input

    Parameters:
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
    return None
