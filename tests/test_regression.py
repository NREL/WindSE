'''
This code will attempt to run all of the regression tests in ./9-Regression/ folder
'''

import pathlib
import pytest
import windse_driver
import os
import yaml
import warnings

### Located Demos ###
home_path = os.getcwd()
reg_path = "../9-Regression/"

### Get Yaml Files ###
yaml_files = sorted(pathlib.Path(__file__, reg_path).resolve().glob('*.yaml'))

### Import the tolerances ###
tolerances = yaml.load(open("tests/9-Regression/Truth_Data/tolerances.yaml"),Loader=yaml.SafeLoader)

###############################################################
######################### Define Tests ########################
###############################################################

### Run Demo Yaml Files
@pytest.mark.parametrize('yaml_file', yaml_files, ids=lambda yaml_file: yaml_file.parts[-2]+"/"+yaml_file.parts[-1])
def test_yaml_execution(yaml_file):

    ### Filter out some benign numpy warnings ###
    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

    ### Run The Windse Simulation
    folder = os.path.split(yaml_file.as_posix())[0]
    os.chdir(folder)
    windse_driver.driver.run_action(params_loc=yaml_file.as_posix())
    os.chdir(home_path)

    ### Grab the name of the run ###
    _, yaml_name = os.path.split(yaml_file)
    yaml_name = yaml_name.split(".")[0]

    ### Import the Truth ###
    truth_loc = folder + "/Truth_Data/" + yaml_name + "_truth.yaml"
    sim_truth = yaml.load(open(truth_loc),Loader=yaml.SafeLoader)

    ### Import the Results ###
    results_loc = folder + "/output/" + yaml_name + "/tagged_output.yaml"
    sim_results = yaml.load(open(results_loc),Loader=yaml.SafeLoader)

    errors = ""

    ### Iterate over the truth and check with the results
    for module_name, truth_dict in sim_truth.items():
        check_dict = sim_results.get(module_name, None)
        tol_dict   = tolerances.get(module_name, {})

        ### Send warning if a module was not checked ###
        if check_dict is None:
            errors += f"Missing Group - {module_name}\n"

        ### Check each value in the module
        for key, truth_value in truth_dict.items():
            check_value = check_dict.get(key,None)

            ### Get test parameters
            tol_value = tol_dict.get(key,[0.0,"absolute"])
            tol = float(tol_value[0])
            check_type = tol_value[1]

            ### Calculate errors ###
            rel_error = abs(check_value-truth_value)/truth_value
            abs_error = abs(check_value-truth_value)

            if check_value is None:
                errors += f"Missing Key - {module_name}: {key} \n"
            
            elif check_type == "absolute" and abs_error > tol:
                errors += f"Value Error - {module_name}: {key} (abs error: {abs_error}, tol: {tol} truth: {truth_value}, check: {check_value})\n"

            elif check_type == "relative" and rel_error > tol:
                errors += f"Value Error - {module_name}: {key} (rel error: {rel_error}, tol: {tol}, truth: {truth_value}, check: {check_value})\n"

    if len(errors)>0:
        errors = yaml_name + "\n" + errors
        raise ValueError(errors)




