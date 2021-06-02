### These must be imported ###
from dolfin import *
from dolfin_adjoint import *

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
    u_ref = solver.problem.bd.bc_velocity
    u     = solver.problem.u_k
    ud = u_ref(x0)-u(x0)
    J = ud[0]
    
    return J
