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
name = "plane_blockage"


### Set default keyword argument values ###
keyword_defaults = {
    "axis": 2,
    "thickness": "rmax",
    "center" : 250,
    "offset_by_mean": False
}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This is a simple blockage metric that integrates the velocity deficit in 
    a plane in front of or above the farm.

    Keyword arguments:
        axis:      the orientation of the plane, "z" for above, "x" for in front
        thickness: how thick of a plane to integrate over
        center:    distance along the axis where the plane is centered 
    '''

    ### Extract keyword arguments
    axis =      int(kwargs.pop("axis"))
    thickness = kwargs.pop("thickness")
    center =    kwargs.pop("center")

    if thickness == "rmax":
        thickness = solver.problem.dom.mesh.rmax()

    ### Get bounds of integration ###
    lb = center - thickness/2.0
    ub = center + thickness/2.0

    # ### Create the Mesh Function to hold the region of integration
    # region = CompiledSubDomain("x[axis]>=lb && x[axis]<=ub", lb=lb, ub=ub, axis=axis)
    # plane_marker = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim())
    # plane_marker.set_all(0)
    # region.mark(plane_marker,1)
    # File("test_"+repr(center)+".pvd")<<plane_marker

    # ### Create measure
    plane_marker = Expression('x[axis] < lb ? 0.0 : (x[axis] > ub ? 0.0 : 1.0)', lb=lb, ub=ub, axis=axis, degree=1)
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
