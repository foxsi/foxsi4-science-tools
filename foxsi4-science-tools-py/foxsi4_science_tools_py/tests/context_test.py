import foxsi4_science_tools_py as f4st

def check_context_load():
    """Test for the context file load in Python. """
    
    # check information to make sure things are loaded properly
    f4_lt = f4st.contextInfo["flight"]["foxsi4"]["launch_time"]["utc"]
    assert f4_lt=="2024-04-17T22:13:00", "Launch time UTC is wrong."

    f4_ltc = f4st.contextInfo["flight"]["foxsi4"]["launch_time"]["clock"]
    assert f4_ltc==0, "Launch clock time is wrong."

if __name__=="__main__":
    check_context_load()