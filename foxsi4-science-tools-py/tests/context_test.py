import yaml
import pathlib

def check_yaml_load():
    """Test for the YAML load in Python. """

    with open(f"{pathlib.Path(__file__).parent}/../../context-information/context.yaml", "r") as file:
        context_info = yaml.safe_load(file)
    
    # check information to make sure things are loaded properly
    f4_lt = context_info["flight"]["foxsi4"]["launch_time"]["utc"]
    assert f4_lt=="2024-04-17T22:13:00", "Launch time UTC is wrong."

    f4_ltc = context_info["flight"]["foxsi4"]["launch_time"]["clock"]
    assert f4_ltc==0, "Launch clock time is wrong."

if __name__=="__main__":
    check_yaml_load()