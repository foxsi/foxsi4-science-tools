# `foxsi4-science-tools` for FOXSI-4  <span>&#129418;</span>

A repository to store all software tools for the FOXSI-4 science.

This will aid in the download of data from co-observing instruments, for example, as well as there analysis. Other code may also be appropriate for here.

**See below on contributing to the repository.**

**Note:** FOXSI-1, -2, and -3 used a completely different system for completely different types of observations compared to FOXSI-4. Therefore, this repository is only appropriate for FOXSI-4.

## Observational Information and the YAML File

<span>&#x1f6a7;</span> The information here is still under construction and may be edited in the future. <span>&#x1f6a7;</span>

Launch Date was 2024 April 17.

| Event               | Date/Time  (UTC) | Time  (AKDT) | Clock Time |
|:--                  |:--               |:-            |:-          |
| Launch Time         | 22:13:00 UTC     | 14:13:00     | T-0        |
| Shutter Door Open   | 22:14:13 UTC     | 14:14:13     | T+73       |
| Observation Start   | 22:14:40 UTC     | 14:14:40     | T+99.647   |
| Shutter Door Closed | 22:20:41 UTC     | 14:20:41     | T+461      |

Times for the flight may only be accurate to within a second. Hi-C launched 1 minute later at 22:14:00 UTC.

YouTube [link](https://www.youtube.com/watch?v=PYM2bRn-5ZY) to countdown audio at Poker Flat Research Range (FOXSI-PI call for count pick-up at timestamp 3:37:42, count picked up at timestamp 3:38:34).

The M1 flare which the sounding rockets were launched on peaked at approximately 2024-04-17T22:08:00 (UTC) and so, for the full flare, the nominal time to investigate is between 2024-04-17T21:57:00$-$22:30:00.

### YAML File

As the repository developes, more information will be added to the [observation parameters YAML](observational-information/observation-parameters.yaml) file.

YAML files are very easy to use in multiple coding languages.

#### In Python

The file's information is readily available in the Python code via

```python
import foxsi4_science_tools_py as f4st

print(f4st.obsInfo)
```

The file can be read in (almost) manually using tools in this repository with

```python
from foxsi4_science_tools_py.io.load_yaml import load_obs_info

obsInfo = load_obs_info()
```

or, [more generally](https://pyyaml.org/wiki/PyYAMLDocumentation),

```python
import yaml

filename = "/path/to/yaml/file/observational-information/observation-parameters.yaml"

with open(filename, "r") as file:
    obsInfo = yaml.safe_load(file)
```

More information on how to access this information in the Python code contained in this repository, see the [README file for the Python package](foxsi4-science-tools-py/foxsi4_science_tools_py/README.md).

#### In IDL

The [syntax](https://www.nv5geospatialsoftware.com/docs/yaml_parse.html) appears to be very easy but seems to require IDLv9 so I have not tested this personally.

```idl
filename = "/path/to/yaml/file/observational-information/observation-parameters.yaml"

obsInfo = yaml_parse(filename)
```

#### In C++

Tools exist, such as `yaml-cpp` found [here](https://github.com/jbeder/yaml-cpp).

## Repository Aim

This repository will hopefully contain all code that proves useful in the working/analysis of general FOXSI-4 science. Several languages might be used and so the top level will be to contain language specific packages.

For example, the `foxsi4-science-tools-py` folder will be a Python package containing all the necessary `.py` files. If other languages are to be included, like `IDL`, then a folder called `foxsi4-science-tools-idl` should be created and used to contain all the `IDL`-ness of the repository. This standard can be applied to the inclusion of other languages (e.g., C++ as `foxsi4-science-tools-cpp`, Shell code as `foxsi4-science-tools-shell`, etc.).

Every `foxsi4-science-tools-<?>` folder should have an "examples" and a "tests" folder. The "examples" folder is a great place to include scripts that show how some of the code in the repository works and the "tests" folder is a fantastic place to put code that tests the tools that have been created.

Additionally, there is also an "examples" and "tests" folder in the top level of the repository so there is a place for anything that fits these folders that spans across code from multiple languages.

## Contributing to the repository

Thank you so much for considering to contribute to the repository! <span>&#127881;</span>

In order to contribute, we ask that you first create your own fork of the repository and then clone that fork to your local machine. Branches of your new fork can be created to develop new features or fix bugs (exciting!). When you are happy with the code in that new branch, a pull request (PR) can be opened which aims to merge the code in your fork's branch into the `main` `foxsi/foxsi4-science-tools` repository. A lot of discussion can be facilitated in an open PR.

**Note:** We aim to _never_ `push` from a local machine to this repository directly. If this happens then it can be very difficult for other contributers to understand what changes are being made and how it affects their own PRs. _If the repository is pushed to directly, in order to help track changes and make them visible to other contributers, the repository will be reverted back to it's state before the push and the undone changes will be asked to be proposed via a PR to then be merged._

## The `external` directory

The `external` directory is a place for software external to the repository. There is no guarentee that it will follow any specific coding language and so it would perhaps not be ideal to place it in a specific coding language `tools` directory. Some may be [`git` submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

To clone submodules when cloning the main repository, try the following:

- `git clone --recurse-submodules https://github.com/foxsi/foxsi4-science-tools.git`
