# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.4.0]

### Added

- Ripp class partially implemented.

### Changed

- Changed structure of main() to accommodate optional input data
- Introduced warning when a RefSeq entry is specified
- Added support for SMILES, compound synonym names, compound evidence, biological activities, biological targets
- Specified allowed entries/regular expressions to Base() class for inheritance
- Moved constants to Base class -> class attributes
- Enabled multiple references per MIBiG entry
- Removed support for patent references



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
