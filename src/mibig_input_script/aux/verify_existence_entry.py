#!/usr/bin/env python3

from pathlib import Path
from typing import Optional, List


def verify_existence_entry(existing: Optional[List]) -> Optional[Path]:
    """Test if MIBiG entry exists in folders and raise error if not

    Parameters:
        existing : `str` OR `None`

    Returns:
        None or Path to entry file
    """
    if existing is not None:
        existing = existing[0]

        mibig_curr_dir_path = (
            Path(__file__).resolve().parent.parent.joinpath("mibig_curr_ver")
        )
        mibig_next_dir_path = (
            Path(__file__).resolve().parent.parent.joinpath("mibig_next_ver")
        )

        if not (Path.exists(mibig_curr_dir_path) and Path.exists(mibig_next_dir_path)):
            raise FileNotFoundError(
                f"""
                One or more of the directories
                '{mibig_curr_dir_path.name}' or '{mibig_next_dir_path.name}'
                containing the MIBiG entries
                does not exist. Check and try running the script again.
                """
            )
        else:
            pass

        if Path.exists(mibig_curr_dir_path.joinpath(existing).with_suffix(".json")):
            return mibig_curr_dir_path.joinpath(existing).with_suffix(".json")
        elif Path.exists(mibig_curr_dir_path.joinpath(existing)):
            return mibig_curr_dir_path.joinpath(existing)
        elif Path.exists(mibig_next_dir_path.joinpath(existing).with_suffix(".json")):
            return mibig_next_dir_path.joinpath(existing).with_suffix(".json")
        elif Path.exists(mibig_next_dir_path.joinpath(existing)):
            return mibig_next_dir_path.joinpath(existing)
        else:
            raise FileNotFoundError(
                f"""
                Could not find the entry '{existing}'.
                Are you sure this entry exists?
                Please check the MIBiG ID for errors and try again.
                """
            )
    else:
        return
