#!/usr/bin/env python3

from importlib import metadata
from mibig_input_script.aux.parse_arguments import parse_arguments
from mibig_input_script.aux.verify_existence_entry import verify_existence_entry

VERSION = metadata.version("mibig-input-script")


def main() -> None:
    """Entry point of the program.

    Serves as main entry point for the program execution.
    It performs the following tasks:
    - Create the program interface via `argparse`
    - Perform checks on
    - Initialize a `MibigEntry()` instance
    - (Optional) Parse an existing MIBiG entry
    - Accept input to add to entry, run verification checks
    - Write the changelog
    - Dump the entry as json file

    Parameters:
        None

    Returns:
        None
    """
    args = parse_arguments(VERSION)

    path_existing = verify_existence_entry(args.existing)

    print(path_existing)


if __name__ == "__main__":
    main()
