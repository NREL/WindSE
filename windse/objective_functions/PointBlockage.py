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
    a single location in the farm.  It does this by creating a small, spherical
    Gaussian centered at the point x0.  The gaussian has a characteristic width
    of mesh.hmax() to ensure it captures enough points for a meaningful 
    measurement of the streamwise velocity.

    Keyword arguments:
        location: where the deficit is evaluated
    '''

    # Get the location of the point measurement and symbolic coordinates for the mesh
    x0 = np.array(kwargs.pop("location"))
    x = SpatialCoordinate(solver.problem.dom.mesh)

    # Calculate dx, dy, dz
    delta_x = x[0] - x0[0]
    delta_y = x[1] - x0[1]
    delta_z = x[2] - x0[2]

    # Calculate the distance to the point, scale by hmax
    distance = (delta_x**2 + delta_y**2 + delta_z**2)/solver.problem.dom.mesh.hmax()

    # Use the distance to generate a sphereical Gaussian function
    spherical_gaussian = exp(-pow(distance, 6.0))

    # Calculate the volume for normalization (this should result in a valid m/s measurement)
    volume = assemble(spherical_gaussian*dx)

    # Use the sphereical Gaussian to measure the streamwise velocity
    if volume <= 1e-10:
        J = np.nan
        print("Warning: No area of integration for point blockage, refine mesh or increase Gaussian width.")
    else:

        ### Evaluate objective ###
        J = assemble(solver.problem.up_k[0]*spherical_gaussian/volume*dx)
    return J
