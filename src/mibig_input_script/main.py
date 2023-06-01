#!/usr/bin/env python3

from __version__ import __version__
from mibig_input_script.aux.parse_arguments import parse_arguments


def main():
    """
    Place to go through the different routes based on args
    - get new information via input
    - load data from existing jsons
    - test for validity
    - dump json files
    """
    parser = parse_arguments(__version__)
    args = parser.parse_args()

    print(args)


if __name__ == "__main__":
    main()
