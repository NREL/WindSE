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
name = "velocity_function"

### Run this code to see which names are already taken ###
# import windse.objective_functions as obj_funcs
# print(obj_funcs.objective_functions.keys())
# print(obj_funcs.objective_kwargs)
# exit()

### Set default keyword argument values ###
# These must be a dictionary and will be passed in via the kwargs.
# Leave empty if no argument are needed. 
keyword_defaults = {
    "velocity_path": None,
}


### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    Description of objective function

    Keyword arguments:
        velocity_path: location of the velocity profile csv file, columns: vel_x, vel_y, vel_z, x, y, z, and a 1 row header
    '''

    # check if the velocity profile is already loaded
    if not hasattr(solver, "velocity_profile_target"):

        # you can extract any kwarg using
        velocity_path = kwargs.pop("velocity_path")
        solver.fprint(f"Loading Velocity: {velocity_path}")

        # Get the function space ###
        V = solver.problem.fs.V

        # Load Velocity Data
        raw_data = np.loadtxt(velocity_path, skiprows=1, delimiter=",").T

        # Create velocity target
        vel = Function(V)

        # reformat velocity data
        stacked_data = np.zeros(vel.vector().get_local().shape)
        stacked_data[0::3] = raw_data[0]
        stacked_data[1::3] = raw_data[1]
        stacked_data[2::3] = raw_data[2]

        # convert velocity data to a dolfin function
        d2v = dof_to_vertex_map(V)
        vel.vector()[:] = stacked_data[d2v]

        # save velocity to the solver
        solver.velocity_profile_target = vel

        # save velocity file 
        solver.params.Save(solver.velocity_profile_target,"velocity_profile_target",subfolder="functions/")


    # Compute difference in velocity
    vel_diff = solver.velocity_profile_target - solver.problem.u_k
    J = assemble(dot(vel_diff, vel_diff)*dx)/solver.problem.dom.volume

    solver.fprint(f"Velocity Difference: {float(J)}")
    # solver.fprint(f"Velocity errornorm: {errornorm(solver.velocity_profile_target,solver.problem.u_k)}")

    return J
