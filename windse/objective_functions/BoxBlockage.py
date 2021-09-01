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
name = "box_blockage"


### Set default keyword argument values ###
keyword_defaults = {
    "p0":   None,
    "p1" :  None,
}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This is a simple blockage metric that integrates the velocity deficit in 
    a plane in front of or above the farm.

    Keyword arguments:
        p0: the corner of the box with the smallest values
        p1: the corner of the box with the largest values
    '''

    ### Extract keyword arguments
    p0 = list(kwargs.pop("p0"))
    p1 = list(kwargs.pop("p1"))

    # ### Create measure
    plane_marker = Expression('((x[0] > xa && x[0] <= xb) && (x[1] > ya && x[1] <= yb) && (x[2] > za && x[2] <= zb)) ? 1 : 0 ', xa=p0[0],ya=p0[1],za=p0[2],xb=p1[0],yb=p1[1],zb=p1[2], degree=2)
    # plane_marker = Expression('x[0] < lb ? 0.0 : (x[0] > ub ? 0.0 : 1.0)', lb=-300.0, ub=300.0, degree=1)
    # plane_marker = Expression('x[0] < pl[0] ? 0.0 : 1.0 ', pl=p0, degree=2)

    # Save the measure (for debugging)
    # test = project(plane_marker,solver.problem.fs.Q)
    # File("test.pvd")<<test
    # exit()

    ### Calculate the volume
    dx = Measure('dx', domain=solver.problem.dom.mesh)
    V = assemble(plane_marker*dx)

    if V <= 1e-10:
        J = np.nan
        print("Warning: No area of integration for plane blockage, refine mesh or increase thickness.")
    else:
        ### Compute velocity deficit
        # u_ref = solver.problem.bd.bc_velocity[0]
        # u     = solver.problem.u_k[0]
        # ud    = (u - u_ref)/u_ref
        # ud    = u

        ### Evaluate objective ###
        J = assemble(plane_marker*solver.problem.u_k[0]*dx)/V
    return J
