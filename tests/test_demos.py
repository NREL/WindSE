import pathlib
import runpy
import pytest
import windse_driver
import glob, os, sys, numpy


### Located Demos ###
home_path = os.getcwd()
demo_path = "../../demo/documented/"

### Get Yaml Files ###
yamls = pathlib.Path(__file__, demo_path+"Yaml_Examples").resolve().glob('*.yaml')

### Get Python Files ###
scripts = pathlib.Path(__file__, demo_path+"Driver_Example").resolve().glob('*.py')

###############################################################
######################### Define Tests ########################
###############################################################

### Run Demo Yaml Files
@pytest.mark.parametrize('yaml', yamls, ids=lambda yaml: yaml.parts[-2]+"/"+yaml.parts[-1])
def test_yaml_execution(yaml):
    folder = os.path.split(yaml.as_posix())[0]
    os.chdir(folder)
    windse_driver.driver.run_action(params_loc=yaml.as_posix())
    os.chdir(home_path)


### Run Demo Python Scripts
@pytest.mark.parametrize('script', scripts, ids=lambda script: script.parts[-2]+"/"+script.parts[-1])
def test_script_execution(script):
    folder = os.path.split(script.as_posix())[0]
    os.chdir(folder)
    runpy.run_path(script.as_posix())
    os.chdir(home_path)