'''
This code will attempt to run all of the regression tests in ./9-Regression/ folder
'''

import pathlib
import pytest
import windse_driver
import os


### Located Demos ###
home_path = os.getcwd()
reg_path = "../9-Regression/"

### Get Yaml Files ###
yamls = sorted(pathlib.Path(__file__, reg_path).resolve().glob('*.yaml'))

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
