This repository contains scripts to help in creating new entries and
modify existing entries for the upcoming version of MIBiG (4.0).

The script takes input, tests it for correctness where possible,
validates the resulting JSON file against the MIBiG JSON schema, and
saves the data as JSON file.

This script is currently for internal use in the Medema/Weber groups
only and requires some knowledge regarding MIBiG entries and command-line usage.

Installation
============

- Create a new virtual env using `python version 3.11`
- Install the package and download requirements using `pip install -e .`

Background
==========

The program is mainly conceptualized to quickly create and modify MIBiG
entries for the next iteration of MIBiG. The script can also be used
for the modification of existing entries (from previous MIBiG iterations).

Currently, the script is limited to create a **minimal** MIBIG entry (or
manipulate the minimum information for a MIBiG entry). This includes
the biosynthetic class, NCBI accession number, organism etc. Support for
additional information is planned for future releases

**Currently** included modules are:
- **base** (essential information)
- **changelog** (writes changelog)
- **write_mibig** (validates and writes MIBiG json files)

**Planned** modules are:
- **compound** (detailed information on compound(s), such as SMILES)
- **ripp** (detailed information on RiPPs)
- **gene_annotation** (annotation of individual genes)


Usage
=====

- Start the program with `python ./src/mibig-input-script/main.py`
- Display the help with the `-h` flag.
- To create a new entry, provide the curator ID together with the `-n` flag (e.g. `... MMZ -n`)
- To modify an existing entry, provide the curator ID, the `-e` flag, and the MIBiG ID (e.g. `... MMZ -e BGC0001234`)

Note 1: The curator ID is taken from `src/mibig_input_script/curators.csv`
Note 2: Start/end coordinates for the locus are mandatory. If the locus equals the BGC, the coordinates still have to be entered ("1" to however long the BGC is).

For developers
==============

In this project, a number of tools are used to keep code and style consistent.
These tools include:
- `black` (v23.3.0)
- `flake8` (v6.0.0)

We recommend using the package `pre-commit` to run these tools before committing.
`pre-commit` can be installed with `pre-commit install`.

Besides, we use type hinting and document code using Google-style docstrings.
A convenient tool to check documentation style is `pycodestyle`.

We use [Semantic Versioning](http://semver.org/) for versioning.

About
=====

## Dependencies

A list of dependencies can be found in [requirements.txt](requirements.txt).


## License

TBA

Authors
=======
- Mitja M. Zdouc (Wageningen University)

