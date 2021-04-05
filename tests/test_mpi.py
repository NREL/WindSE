'''
This code will attempt to run all of the regression tests in ./9-Regression/ folder
'''

import pathlib
import pytest
import os, sys
import yaml
import warnings
import subprocess

### Located Demos ###
home_path = os.getcwd()
reg_path = "../9-Regression/"

### Get Yaml Files ###
yaml_files = sorted(pathlib.Path(__file__, reg_path).resolve().glob('*_Unsteady.yaml'))

### Import the tolerances ###
tolerances = yaml.load(open("tests/9-Regression/Truth_Data/tolerances.yaml"),Loader=yaml.SafeLoader)

### Get current status of modules
default_modules = sys.modules.keys()

## Set the number of processors to test with
num_procs = 2

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
    # from windse_driver import driver

    ### Grab the name of the run ###
    _, yaml_name = os.path.split(yaml_file)
    yaml_name = yaml_name.split(".")[0]
    parallel_yaml_name = yaml_name + "_np_%d" % (num_procs)

    ### Set the number of processes and launch using mpirun as if from the command line
    cl_command = "mpirun -n %d windse run %s -p general:name:%s" % (num_procs, yaml_file, parallel_yaml_name)
    # print(cl_command)
    cl_output = subprocess.run(cl_command.split())

    # driver.run_action(params_loc=yaml_file.as_posix())
    os.chdir(home_path)

    ### Import the Truth ###
    truth_loc = folder + "/Truth_Data/" + yaml_name + "_truth.yaml"
    sim_truth = yaml.load(open(truth_loc),Loader=yaml.SafeLoader)

    ### Import the Results ###
    results_loc = folder + "/output/" + parallel_yaml_name + "/tagged_output.yaml"
    sim_results = yaml.load(open(results_loc),Loader=yaml.SafeLoader)

    errors = ""

    ### Iterate over the truth and check with the results
    for module_name, truth_dict in sim_truth.items():
        check_dict = sim_results.get(module_name, None)
        tol_dict   = tolerances.get(module_name, {})

        ### Send warning if a module was not checked ###
        if check_dict is None:
            errors += f"Missing Group - {module_name}\n"

        ### Get Wildcard Tolerances for module ###
        wild_tol_keys = []
        for key in tol_dict.keys():
            if "*" in key:
                wild_tol_keys.append(key)

        ### Check each value in the module
        for key, truth_value in truth_dict.items():
            check_value = check_dict.get(key,None)
            tol_value = None

            ### Check if there is a valid wildcard ###
            use_wildcard = False
            for wild_key in wild_tol_keys:
                filtered_wild_key = wild_key.replace('*', '')
                if filtered_wild_key in key:
                    wild_tol_value = tol_dict[wild_key]
                    use_wildcard = True

            ### Check if the exact key is available ###
            if key in tol_dict.keys():
                tol_value = tol_dict[key]

            ### Check if wildcard tolerance key is available ###
            elif use_wildcard:
                tol_value = wild_tol_value

            ### Set the default tolerances for a float ###
            elif isinstance(check_value,float):
                tol_value = [1e-4,"absolute"]

            ### Set the default tolerances for a float ###
            elif isinstance(check_value,int):
                tol_value = [0,"absolute"]

            ### Get tolerance parameters ###
            tol = float(tol_value[0])
            check_type = tol_value[1]


            if check_value is None:
                errors += f"Missing Key - {module_name}: {key} \n"
            else:
                ### Calculate errors ###
                abs_error = abs(check_value-truth_value)
                if truth_value != 0:
                    rel_error = abs(check_value-truth_value)/truth_value

            if check_type == "absolute" and abs_error > tol:
                errors += f"Value Error - {module_name}: {key} (abs error: {abs_error}, tol: {tol} truth: {truth_value}, check: {check_value})\n"

            elif check_type == "relative" and rel_error > tol:
                errors += f"Value Error - {module_name}: {key} (rel error: {rel_error}, tol: {tol}, truth: {truth_value}, check: {check_value})\n"

    if len(errors)>0:
        errors = parallel_yaml_name + "\n" + errors
        raise ValueError(errors)




