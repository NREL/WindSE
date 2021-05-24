#######################################################################
####################### Preamble (do not edit) ########################
#######################################################################

import __main__
import os

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

#######################################################################
#######################################################################

# Feel free to rename this file anything. if the file starts with a '_',
# it will be ignored in the import.

### Additional import statements ###
import numpy as np
import math

### Declare Unique name
name = "objective_name"

### Run this code to see which names are already taken ###
# import windse.objective_functions as obj_funcs
# print(obj_funcs.objectives_dict.keys())
# exit()

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False):
    '''
    This is the function that you use to calculate the objective.
    Do not rename this function. The solver input should contain 
    everything you need. For example: solver.problem.u_k is the 
    most recent velocity field. You can always check what 
    attributes an object has by running print(dir(object)). For 
    example print(dir(solver.problem))) will display every 
    attribute solver.problem has. Avoid adding attributes to 
    solver unless you really, REALLY know what your are doing.

    The inflow_angle input is the angle off the x-axis (different 
    from yaw). This is a legacy option and probably should not be
    used as it will likely be removed in the future.

    Finally, first_call is a flag that let this function know if 
    this is the first time it has been called. This is useful for
    saving objective specific files as it lets you know you need
    to create the file if true and append to the file if false. 

    To make sure your function plays nice with dolfin_adjoint, 
    try to exclusively use dolfin operation such as assemble and 
    dolfin's versions of sqrt, sin(), cos(), etc. Avoid using 
    external operations from numpy or other packages when dealing 
    with any value/function directly tied to  the control 
    variables. In the end check the type of J. With  
    dolfin_adjoint imported, type(J) should return something like
    pyadjoint.adjfloat.AdjFloat or possibly something from 
    numpy_adjoint. Either way, it should probably reference 
    'adjoint' Note: to check this, make sure 
    general:dolfin_adjoint is set to True in the yaml file.
    '''

    ### this is a dummy value and can be removed
    J = Constant(0.0)

    return J
