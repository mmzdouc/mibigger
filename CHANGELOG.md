# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.5.4] - 03-10-2023

### Changed

- Made "completeness" a fields that needs to be defined explicitly.

## [0.5.1] - 10-09-2023

### Changed

- Fixed wrongly formatted SMILES (double-escaped characters)

## [0.5.0] - 10-09-2023

### Added

- Pre-commit now includes the `check_valid.py` script to test json files

### Changed

- Renamed the project to `mibigger`
- Updated `schema.json` to version 2.13
- Renamed folder `aux` to `helper` (`aux` is a forbidden dir name in Windows)

## [0.4.1] - 23-06-2023

### Added

- Added self-updating options for tailoring reactions and ripp subclasses

### Changed

- Added additional evidence options to "loci"; modified JSON schema file


## [0.4.0] - 20-06-2023

### Added

- Added optional input functions for Gene Annotation.
- Added optional input functions for RiPPs.
- Added optional input functions for compounds.

### Changed

- Changed structure of main() to accommodate optional input data
- Introduced warning when a RefSeq entry is specified
- Specified allowed entries/regular expressions to Base() class for inheritance
- Moved constants to Base class -> class attributes
- Enabled multiple references per MIBiG entry
- Removed support for patent references as publication


## [0.3.0] 12-06-2023

### Changed

- Renamed class Base to MibigEntry.
- The MibigEntry class now generates a dict called mibig_entry that is consecutively modified by other classes.
- Class Changelog now modifies the mibig_entry dict instead of generating a new dict.
- Class WriteMibig now uses the mibig_entry dict directly instead of concatenating
    the dicts into a new dict.

## [0.2.0] 11-06-2023

### Added

- Added duplicate testing to class WriteMibig.

## [0.1.0] 09-06-2023

### Added

- Initial version containing the base, changelog, and export modules.
- Automated validation using the MIBiG JSON schema.
