# `foxsi4-science-tools-py`

A repository to store all software tools for FOXSI-4 science in Python (woo, Python).

The tools being developed in Python should be placed in the `foxsi4_science_tools_py` folder (note the underscores and not dashes).

There is an "examples" and a "tests" folder. The "examples" folder is a great place to include scripts that show how some of the code in the repository works and the "tests" folder is a fantastic place to put code that tests the tools that have been created.

More information will be placed here with regards as to how this package is recommended to be used.

## Install tips

In order to work with some preliminary data, it would be instructive to set up a virtual environment (more information below) and install some packages needed. One way to do this is to:

1. Create an environment with `conda create -n foxsi4-science-tools-env python`.
2. After activating that environment, and in the directory with the `setup.py` file, the command `pip install -e .` should be possible to install the Python code into that environment.

## A useful Python tip

It might be a good idea to look into ([conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)) virtual environments if you are not familiar, this includes looking into them yourself or getting in touch with someone to help explain.
