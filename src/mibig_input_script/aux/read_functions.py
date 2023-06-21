"""
Collection of auxiliary read functions.

Contains auxiliary functions that read from files.
"""
import csv
import json
from pathlib import Path
import pandas as pd
from typing import Dict, List


def read_mibig_json(path: Path | None) -> Dict | None:
    """Read an existing MIBiG json file.

    Parameters:
        `path` : a `Path` object OR `None`

    Returns:
        MIBiG entry in form of `dict` OR `None`
    """
    if path is not None:
        with open(path.resolve(), "r") as jsonfile:
            jsonfile_contents = jsonfile.read()
            return json.loads(jsonfile_contents)
    else:
        return None


def get_curators(ROOT: Path) -> Dict:
    """Get a Dict of possible curators and their IDs defined in `../curators.csv`.

    Parameters:
        `ROOT`: `Path` object indicating "root" directory of script

    Returns:
        Returns a dict of curator ids
    """
    csv_path = ROOT / "curators.csv"
    df = pd.read_csv(csv_path)
    return dict(zip(list(df.loc[:, "initials"]), list(df.loc[:, "id"])))


def get_curator_email(ROOT: Path, curator: str) -> str:
    """Get the email of the specified curator.

    Parameters:
        `ROOT`: `Path` object indicating "root" directory of script
        `curator` : `str` indicating the initials of the curator

    Returns:
        Returns an email as string
    """
    csv_path = ROOT / "curators.csv"
    df = pd.read_csv(csv_path)
    return df.loc[df["initials"] == curator]["email"].iloc[0]
