import pathlib
import runpy
import pytest
import glob, os, sys, numpy
import copy
from importlib import reload


### Located Demos ###
home_path = os.getcwd()
demo_path = "../../demo/documented/"

### Get Yaml Files ###
yamls = sorted(pathlib.Path(__file__, demo_path+"Yaml_Examples").resolve().glob('*.yaml'))

### Get Python Files ###
drivers = sorted(pathlib.Path(__file__, demo_path+"Driver_Example").resolve().glob('*.py'))

### Get current status of modules
default_modules = list(sys.modules.keys())

###############################################################
######################### Define Tests ########################
###############################################################

### Run Demo Drivers
@pytest.mark.parametrize('driver', drivers, ids=lambda driver: driver.parts[-2]+"/"+driver.parts[-1])
def test_driver_execution(driver):
    folder = os.path.split(driver.as_posix())[0]
    os.chdir(folder)
    runpy.run_path(driver.as_posix())
    os.chdir(home_path)

### Run Demo Yaml Files
@pytest.mark.parametrize('yaml', yamls, ids=lambda yaml: yaml.parts[-2]+"/"+yaml.parts[-1])
def test_yaml_execution(yaml):
    folder = os.path.split(yaml.as_posix())[0]
    os.chdir(folder)
    import windse_driver
    windse_driver.driver.run_action(params_loc=yaml.as_posix())
    os.chdir(home_path)
