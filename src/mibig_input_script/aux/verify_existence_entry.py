from pathlib import Path
from typing import List


def verify_existence_entry(
    existing: List | None,
    ROOT: Path,
) -> Path | None:
    """Test if MIBiG entry exists in folders and raise error if not

    Parameters:
        `existing` : `List` OR `None`
        `ROOT` : `Path` object indicating "root" directory of script

    Returns:
        `Path` OR `None` : a Path object pointing towards an MIBiG json file.
    """
    if existing is not None:
        mibig_entry = existing[0]

        mibig_curr_dir_path = ROOT.joinpath("mibig_curr_ver")
        mibig_next_dir_path = ROOT.joinpath("mibig_next_ver")

        if not (Path.exists(mibig_curr_dir_path) and Path.exists(mibig_next_dir_path)):
            message = (
                f"One or more of the directories\n"
                f"'{mibig_curr_dir_path.name}' or '{mibig_next_dir_path.name}\n"
                f"containing the MIBiG entries do not exist.\n"
                f"Check your installation and try running the script again.\n"
            )
            raise FileNotFoundError(message)
        else:
            pass

        if Path.exists(mibig_curr_dir_path.joinpath(mibig_entry).with_suffix(".json")):
            return mibig_curr_dir_path.joinpath(mibig_entry).with_suffix(".json")
        elif Path.exists(mibig_curr_dir_path.joinpath(mibig_entry)):
            return mibig_curr_dir_path.joinpath(mibig_entry)
        elif Path.exists(
            mibig_next_dir_path.joinpath(mibig_entry).with_suffix(".json")
        ):
            return mibig_next_dir_path.joinpath(mibig_entry).with_suffix(".json")
        elif Path.exists(mibig_next_dir_path.joinpath(mibig_entry)):
            return mibig_next_dir_path.joinpath(mibig_entry)
        else:
            message = (
                f"Could not find the entry:\n"
                f"'{mibig_entry}'\n"
                f"Are you sure this entry exists? Please check again.\n"
                f"If you want to create a new entry, please use the '-n' flag.\n"
                f"New entries automatically receive a new MIBiG accession number.\n"
            )
            raise FileNotFoundError(message)
    else:
        return None
