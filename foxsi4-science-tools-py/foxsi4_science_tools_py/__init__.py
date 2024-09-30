"""
Top level of the `foxsi4-science-tools-py` package.

Namespace
---------

The base `foxsi4_science_tools_py` namespace includes:

`~foxsi4_science_tools_py.obsInfo`
    The information stored in a YAML file that includes information from
    the flight and also the flare that was observed.

Examples
--------
# importing the module is as easy as:
>>> import foxsi4_science_tools_py as f4st

# then accessing, e.g., the observational information for the flight/flare
>>> print(f4st.obsInfo)
{'flight': {'foxsi4': {'launch_time': {'utc': '2024-04-17T22:13:00', 'akdt':...}

# that can be accessed like a native Python dictionary
>>> print(f4st.obsInfo["flight"]["foxsi4"]["launch_time"]["clock"])
0
"""

from foxsi4_science_tools_py.io.load_yaml import load_obs_info

# for global observational info
obsInfo = load_obs_info()