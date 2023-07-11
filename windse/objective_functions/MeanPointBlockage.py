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

### Declare Unique name
name = "mean_point_blockage"


### Set default keyword argument values ###
keyword_defaults = {
    "z_value": 240
}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This is a simple blockage metric that evaluates the velocity deficit at 
    a single location above the mean location of all turbine in the farm.

    Keyword arguments:
        z_value: z location to evaluate
    '''
    farm_x, farm_y, farm_z = solver.problem.farm.get_hub_locations().T
    x0 = np.mean(farm_x)
    y0 = np.mean(farm_y)
    z0 = float(kwargs.pop("z_value"))
    p0 = np.array([x0,y0,z0])

    # u_ref = solver.problem.bd.bc_velocity
    # u     = solver.problem.u_k
    # ud = (u(p0)-u_ref(p0))/u_ref(p0)
    # ud = u(p0)
    J = solver.problem.up_k(p0)[0]
    
    return J
