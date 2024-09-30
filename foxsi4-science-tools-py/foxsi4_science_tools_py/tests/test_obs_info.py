import foxsi4_science_tools_py as f4st

def test_load_obs_info():
    """Test for the observational infor file load in Python. """
    
    # check information to make sure things are loaded properly
    f4_lt = f4st.obsInfo["flight"]["foxsi4"]["launch_time"]["utc"]
    assert f4_lt=="2024-04-17T22:13:00", \
        "Launch time UTC is wrong in global observational information."

    f4_ltc = f4st.obsInfo["flight"]["foxsi4"]["launch_time"]["clock"]
    assert f4_ltc==0, \
        "Launch clock time is wrong in global observational information."

if __name__=="__main__":
    test_load_obs_info()