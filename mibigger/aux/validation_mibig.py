"""Validate MIBiG entry.

Validates MIBig entry by checking against JSON schema and performing
futher checks.
"""
from pathlib import Path
from mibigger.classes.class_write_mibig import WriteMibig

import json
import jsonschema


def validation_mibig(write_mibig: WriteMibig, ROOT: Path) -> None:
    """Validate the MIBiG entry before writing it to disk.

    Parameters:
        `write_mibig`: an instance of the class WriteMibig
        `ROOT`: an instance of the class Path

    Returns:
        None
    """
    json_string = write_mibig.return_json_string()

    with open(ROOT.joinpath("aux").joinpath("schema.json")) as schema_handle:
        schema = json.load(schema_handle)
    try:
        jsonschema.validate(json.loads(json_string), schema)
    except jsonschema.ValidationError as e:
        print("MIBiG JSON is invalid. Error:", e)
        print("Abort file storage.")
        exit()

    return
