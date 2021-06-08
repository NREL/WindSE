### These must be imported ###
from dolfin import *
from dolfin_adjoint import *

### Additional import statements ###
import numpy as np

### Declare Unique name
name = "plane_blockage"


### Set default keyword argument values ###
keyword_defaults = {
    "axis": 2,
    "thickness": 50,
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

    ### Get bounds of integration ###
    lb = center - thickness/2.0
    ub = center + thickness/2.0
    region = CompiledSubDomain("x[axis]>=lb && x[axis]<=ub", lb=lb, ub=ub, axis=axis)

    ### Create the Mesh Function to hold the region of integration
    plane_marker = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim())
    plane_marker.set_all(0)
    region.mark(plane_marker,1)

    ### Create measure
    dx = Measure('dx', subdomain_data=plane_marker)

    ### Compute velocity deficit
    u_ref = solver.problem.bd.bc_velocity[0]
    u     = solver.problem.u_k[0]
    ud    = u_ref - u

    ### Evaluate objective ###
    J = assemble(ud*dx(1))
    return J
