################## HEADER DO NOT EDIT ##################
import os
import __main__

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    from dolfin import *
    from dolfin_adjoint import *
########################################################

# Feel free to rename this file anything. if the file starts with a '_',
# it will be ignored in the import.

### Additional import statements ###
import numpy as np

### Declare Unique name
name = "objective_name"

### Run this code to see which names are already taken ###
# import windse.objective_functions as obj_funcs
# print(obj_funcs.objective_functions.keys())
# print(obj_funcs.objective_kwargs)
# exit()

### Set default keyword argument values ###
# These must be a dictionary and will be passed in via the kwargs.
# Leave empty if no argument are needed. 
keyword_defaults = {}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
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

    The first_call flag lets this function know if this is the 
    first time it has been called. This is useful for saving 
    objective specific files as it lets you know you need to 
    create the file if true and append to the file if false. 

    Any additional keyword argument can be popped from
    kwargs. These arguments are loaded through the yaml file and
    directly passed to the objective() function.

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

    Finally, the output of the objective function should be a 
    single value.

    When creating a new objective function make sure to fill out
    the doc string below and delete everything above this line:
    ------------------------------------------------------------

    Description of objective function

    Keyword arguments:
        kw1: description of argument
        kw2: description of argument
    '''

    ### you can extract any kwarg using ###
    # value = kwargs.pop("value_name")

    ### this is a dummy value and can be removed
    J = Constant(0.0)

    return J
