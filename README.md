This repository contains scripts to help in creating new entries for
the upcoming version of MIBiG (4.0).

This script is currently for internal use only and requires some
knowledge regarding MIBiG entries and command-line usage.

Installation
============

- Create a new virtual env using `python version 3.11`
- Install the requirements using `pip install -e .`
- Set up pre-commit via `pre-commit install`

Usage
=====

- The program can be started with `python ./src/mibig-input-script/main.py`
- The help can be displayed using the `-h` flag (i.e. `python ./src/mibig-input-script/main.py` )
- The program allows to create new entries OR manipulate existing entries. These two actions are mutually exclusive
- For both creating/manipulating an entry, the curator need to "log in" using their curator ID as first command line argument (i.e. `python ./src/mibig-input-script/main.py MMZ`)
- New entries must be initiated using the flags `-n --minimal`. If necessary, additional flags can be set to add more information (e.g. `--gene_ann`
- Existing entries can be manipulated using the `-e <MIBIG_ACCESSION_NUMBER>`, using additional flags specifying the type of information that is going to be manipulated (e.g. `minimal`, `gene_ann`). These additional flags


About
=====

## Dependencies

A list of dependencies can be found in [requirements.txt](requirements.txt).

## License

TBA

For developers
==============

## Contributing

We document code with Google-style docstrings -> see the examples below

```
class MyClass:
    """Brief summary of the class.

    More detailed description of the class.

    Attributes:
        attr1 (type): Description of attr1.
        attr2 (type): Description of attr2.

    Methods:
        method1(arg1 : type, arg2 : type) -> type:
            Description of method1.
        method2(arg1 : type) -> type:
            Description of method2.

    """
    # class implementation goes here
```

```
def myfunct(arg1: float, arg2: float) -> float:
    """Brief summary of the function.

    Args:
        arg1: Description of arg1.
        arg2: Description of arg2.

    Returns:
        Description of return value

    Notes:
        Optional details on function
    """
    # function implementation goes here
```


## Versioning

We use [Semantic Versioning](http://semver.org/) for versioning.

## TO DO (REMOVE LATER)

- Add pre-commit code from Simon to main.py
