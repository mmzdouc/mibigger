Overview
========

This repository contains scripts to help in creating new entries and
modify existing entries for the upcoming version of MIBiG (4.0).

The script takes input, tests it for correctness where possible,
validates the resulting JSON file against the MIBiG JSON schema, and
saves the data as JSON file.

This script is currently for internal use in the Medema/Weber groups
only and requires some knowledge regarding MIBiG entries, use of git, and command-line
usage.

Curation Strategy
=================

Curation is a prime example of distributed development: multiple people work
on the same project, adding new information or editing existing ones. To prevent
conflicting files and parallel edits, the following strategies are implemented:
- Branching model: all curators work on separate branches, which branch out of the
  main (production) branch. Once a (few) new entries are created/modified, the
  curator issues a pull request to the main branch. After successful review and
  addition of changes, the branch is incorporated in the main.
- Reviewing: each pull request needs to be reviewed to guarantee factual correctness
  (four-eyes-principle). During the pull request, a reviewer can be assigned. For a
  list of reviewers and their specialities, see below.
- Communication: version control is no replacement for communication. Both curators
  and reviewers need to communicate which entries they are working on to prevent
  possible redundant work. Please announce your work via the Slack channel.
- Entry naming: To prevent identical names, all new entries are automatically
  prefixed with the curators initials (or MIBiG annotator ID, if already available)
- Entry validation: each entry is validated using the `validated.py` script of Kai
  Blin. The validation can be automated by using the package `pre-commit` (see below).


Download, Installation, Curation
================================

- Clone `mibigger` to create a local version
- Create and activate a new virtual env using `python version 3.11` (e.g. with `conda`)
- Install the package and download requirements using `pip install -e .`
- Install pre-commit using `pre-commit install`
- Create a new branch using `git checkout -b dev_mmz` (with your initials instead
  `mmz`)
- Push your new branch to GitHub using `git push origin dev_mmz`
- Set the remote branch to push to using `git branch --set-upstream-to=origin/dev_mmz
  dev_mmz`
- Add your initials/curator_id/name/email to the list of curators in
  `src/mibig_input_script/curators.csv`
- Add/modify entries (see below) and commit to your branch. Pre-commit will
  automatically check committed .json-files for correctness and prevent incorrect
  entries.
- Ideally, after each edit/new entry, the changes should be pulled into the main.
  The easiest way of doing so is via the GitHub webpage. Create a pull request and
  request reviewer(s) from the list of collaborators.
- The reviewer (one of us) will double-check the entry for consistency  (see
  instructions below).
- Once reviewed and fixed, the new entries are merged into the main branch.

Usage
=====

- Start the program with `python ./mibigger/main.py`
  - Display the help with the `-h` flag.
  - To create a new entry, provide the curator ID together with the  `-n` flag (e.g.
    `python ./mibigger/main.py MMZ -n`)
  - To modify an existing entry, provide the curator ID, the `-e` flag, and the
    MIBiG ID (e.g. `python ./mibigger/main.py MMZ -e BGC0001234`)

Curation Tips
=================

- **Communication**: In distributed development, version control is no substitute
  for communication.  Please briefly announce the entries you are working on in the
  Slack  channel and check if anybody is work on the entry. Do frequent
  pull requests to keep the main up-to-date.
- **Typos**: please double-check your additions for typos to make it  easier for
  the reviewers. Although we perform many checks for erroneous data, we can't catch
  all mistakes.
- **Coordinates**: Start/end coordinates for the locus are mandatory.  If the locus
  equals the BGC, the coordinates still have to be entered  ("1" to "n",  long the
  BGC is).
- **Preprints**: BGCs described in preprints are okay to add but additional  care
  has to be taken regarding the authenticity of data. Please also mention the signal
  word **preprint** in the comment - this way, we can identify the  preprint-based
  entries in the future.
- **GenBank/RefSeq**: GenBank accessions are strongly preferred over RefSeq
  accessions since the latter can be changed abruptly. Sometimes, the GenBank  entry
  cannot be taken due to missing annotations. In these cases, the RefSeq entry is
  tolerated.

Reviewers (TBA)
=========

- Barbara: NRPS, PKS
- David: NRPS, PKS
- Frida: PKS, RiPPs
- Jorge: All fungal-related
- Marnix: All
- Mitja: RiPPs


Reviewer instructions
=====================

The input program performs several checks to validate new entries. Still, mistakes
can happen and the review can help to catch them:

1) Are there any obvious typos or mistakes?
2) Is the NCBI accession correct? Are the coordinates correct? Search them in the
   NCBI (e.g. by replacing `ACCESSION`, `START`, `STOP` in this URL:
   `https://www.ncbi.nlm.nih.gov/nuccore/ACCESSION?from=START&to=STOP`)
3) If the entry is considered okay, "Request changes" and add your curator ID,
   including a brief comment  (e.g. reviewed by XYZ)
4) The person requesting the review then adds the reviewer to the changelog with an
   "empty edit" and pushes the changes to the repo.
5) The reviewer approves and merges the new entry into the main.

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
