import foxsi4_science_tools_py as f4st

def test_load_yaml():
    """Test for the observational info file load in Python. """

    obs_info = f4st.io.load_yaml.load_obs_info()
    
    # check information to make sure things are loaded properly
    f4_lt = obs_info["flight"]["foxsi4"]["launch_time"]["utc"]
    assert f4_lt=="2024-04-17T22:13:00", \
        "Launch time UTC is wrong from loading YAML."

    f4_ltc = obs_info["flight"]["foxsi4"]["launch_time"]["clock"]
    assert f4_ltc==0, \
        "Launch clock time is wrong from loading YAML."

if __name__=="__main__":
    test_load_yaml()