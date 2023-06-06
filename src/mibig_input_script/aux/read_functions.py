#!/usr/bin/env python3

import csv
import json
from pathlib import Path
from typing import Dict, List


def read_mibig_json(path: Path | None) -> Dict | None:
    """Read an existing MIBiG json file

    Parameters:
        `path` : a `Path` object OR `None`

    Returns:
        `json` MIBiG entry in form of `dict` OR `None`
    """
    if path is not None:
        with open(path.resolve(), "r") as jsonfile:
            jsonfile_contents = jsonfile.read()
            return json.loads(jsonfile_contents)
    else:
        return None


def get_curators(ROOT: Path) -> List[str]:
    """Get the list of possible curators defined in ../curators.csv

    Parameters:
        `ROOT`: `Path` object indicating "root" directory of script

    Returns:
        Returns a list of curator ids
    """
    csv_path = ROOT / "curators.csv"

    with open(csv_path, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        return [row[0] for row in csvreader]
