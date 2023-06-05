#!/usr/bin/env python3

import json
from pathlib import Path
from typing import Dict


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
