This repository contains scripts to help in creating new entries and
modify existing entries for the upcoming version of MIBiG (4.0).

The script takes input, tests it for correctness where possible,
validates the resulting JSON file against the MIBiG JSON schema, and
saves the data as JSON file.

This script is currently for internal use in the Medema/Weber groups
only and requires some knowledge regarding MIBiG entries and command-line usage.

Download, Installation, and Curation Strategy
============

- Create a new branch with your initials as suffix (e.g. `dev_mmz`) or fork the repo
- Clone the new branch/fork locally
- Create a new virtual env using `python version 3.11` (e.g. with `conda`)
- Install the package and download requirements using `pip install -e .`
- Add your initials(curator_id/name/email to the list of curators in `src/mibig_input_script/curators.csv`
- Add/modify entries (see below)
- Once you have added some entries, create a pull request into main and request a reviewer
- The reviewer (one of us) will double-check the entry for consistency
- Once reviewed, the new entries are merged into the main branch

Usage
=====

- Start the program with `python ./src/mibig-input-script/main.py`
- Display the help with the `-h` flag.
- To create a new entry, provide the curator ID together with the `-n` flag (e.g. `... MMZ -n`)
- To modify an existing entry, provide the curator ID, the `-e` flag, and the MIBiG ID (e.g. `... MMZ -e BGC0001234`)

Curation Strategy
=================

- **Communication**: In distributed development, version control is no substitute for communication. Please briefly announce the entries you are working on in the Slack channel and check if anybody is work on the entry. Do frequent pushes and pull requests to keep the main up-to-date with your work.
- **Typos**: please double-check the data you add for typos to make it easier for the reviewers. Altough we perform many checks for erroneous data, we can't catch all mistakes.
- **Coordinates**: Start/end coordinates for the locus are mandatory. If the locus equals the BGC, the coordinates still have to be entered ("1" to however long the BGC is).
- **Preprints**: BGCs described in preprints are okay to add but additional care has to be taken regarding the authenticity of data. Please also mention the signal word "preprint" in the comment - this way, we can identify the preprint-based entries in the future.
- **GenBank/RefSeq**: GenBank accessions are strongly preferred over RefSeq accessions since the latter can be changed abruptly. Sometimes, the GenBank entry cannot be taken due to missing annotations. In these cases, the RefSeq entry is tolerated.

Background
==========

The program is mainly conceptualized to quickly create and modify MIBiG
entries for the next iteration of MIBiG. The script can also be used
for the modification of existing entries (from previous MIBiG iterations).
The script is capable to create minimal entries as well as add annotations
on RiPP precursors and gene annotation. Other functionalities, such as
on NRPS or PKS annotation is not included yet.

**Currently** included modules are:
- **Base** (defines constants and allowed entries)
- **Changelog** (writes changelog)
- **MibigEntry** (writes the essential information of a MIBiG entry)
- **Genes** (writes gene annotation data)
- **RiPP** (writes RiPP annotation data)
- **WriteMibig** (validates and writes MIBiG json files)

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

MIT license (see [LICENSE](LICENSE.md))

Authors
=======

- Mitja M. Zdouc (Wageningen University)
