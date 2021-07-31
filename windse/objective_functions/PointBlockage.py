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
name = "point_blockage"


### Set default keyword argument values ###
keyword_defaults = {
    "location": (0,0,0),
    "offset_by_mean": False
}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This is a simple blockage metric that evaluates the velocity deficit at 
    a single location in the farm.

    Keyword arguments:
        location: where the deficit is evaluated
    '''

    x0 = np.array(kwargs.pop("location"))
    # u_ref, v_ref, w_ref = solver.problem.bd.bc_velocity.split()
    # u,     v,     w     = solver.problem.u_k.split()
    # u = solver.problem.p_k
    # ud = (u(x0)-u_ref(x0))/u_ref(x0)
    # ud = u(x0)
    J = solver.problem.up_k(x0)[0]

    return J
