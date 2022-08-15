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

### Additional import statements ###
import numpy as np
import math
import os

### Declare Unique name
name = "power"

### Set default keyword argument values ###
keyword_defaults = {}

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, annotate=True, **kwargs):
    '''
    The "power" objective function calculates the power using actuator disks, 
    by computing the integral of turbine force dotted with velocity. 
    '''

    J = solver.problem.farm.compute_power(solver.problem.u_k,inflow_angle)


    J *= solver.params.hard_scaling_factor

    return J
