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
name = "cyld_kernel"

### Run this code to see which names are already taken ###
# import windse.objective_functions as obj_funcs
# print(obj_funcs.objective_functions.keys())
# print(obj_funcs.objective_kwargs)
# exit()

### Set default keyword argument values ###
# These must be a dictionary and will be passed in via the kwargs.
# Leave empty if no argument are needed. 
keyword_defaults = {'type': 'above', # 'upstream'
                    'radius': 0.5,
                    'length': 3.0,
                    'sharpness': 6,
                    }

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    This is a blockage metric that measures the velocity within a 
    Gaussian cylinder located upstream from each turbine and aligned
    with the rotor's rotational axis or overhead from each turbine
    and aligned with the mast.  The cylindrical Gaussian field is
    formed by intersecting a radial Gaussian and a streamwise Gaussian.

    Keyword arguments:
        type:      The orientation of the Gaussian cylinder, 
                   "upstream" for Gaussians shifted to measure the velocity
                   directly upstream from each turbine, "above" to orient the
                   Gaussians over the top of each turbine.
        radius:    The radius of the cylinder, expressed in units of rotor diameter (RD).
                   The default of 0.5 sets the radius to match the turbine radius, 0.5*RD
        length:    The length of the cylinder, expressed in units of rotor diameter (RD).
                   The default of 3.0 means the measurement volume has an axial length of 3*RD.
        sharpness: The sharpness value controls the severity with which the Gaussian field 
                   drops to zero outside the volume of the Gaussian cylinder. The sharpness
                   value *must* be an even number.  Smaller values result in smoother transitions
                   over longer length scales, larger values result in more abrupt transitions.
                   For very large values, the transition becomes a near step change which may
                   not have valid values of the derivative.  The default setting of 6 is 
                   a good starting point.
    '''

    def rotate_and_shift_points(x, x0, yaw, kp, RD):

        xs =  cos(yaw)*(x[0]-x0[0]) + sin(yaw)*(x[1]-x0[1])
        ys = -sin(yaw)*(x[0]-x0[0]) + cos(yaw)*(x[1]-x0[1])

        if x.geometric_dimension() == 3:
            zs = x[2] - x0[2]
        else:
            zs = 0.0

        if kp['type'] == 'upstream':
            xs = xs + 0.5*kp['length']
        elif kp['type'] == 'above':
            zs = zs - RD - 0.5*kp['length']

        return [xs, ys, zs]

    def build_cylindrical_kernels(solver, kp):
        x = SpatialCoordinate(solver.problem.dom.mesh)

        kernel_exp = 0

        kernel_exp_list = []

        for k in range(solver.problem.farm.numturbs):

            # Get a convenience copy of turbine k's location
            mx = solver.problem.farm.turbines[k].mx
            my = solver.problem.farm.turbines[k].my
            mz = solver.problem.farm.turbines[k].mz
            RD = solver.problem.farm.turbines[k].RD
            x0 = [mx, my, mz]

            # Get a convenience copy of turbine k's yaw
            yaw = solver.problem.farm.turbines[k].myaw

            xs = rotate_and_shift_points(x, x0, yaw, kp, RD)

            if kp['type'] == 'upstream':
                # Place the cylinders upstream from the rotor aligned with the hub axis
                ax = xs[0]/(0.5*kp['length'])
                axial_gaussian = exp(-pow(ax, kp['sharpness']))
                rad = (xs[1]**2 + xs[2]**2)/kp['radius']**2
                radial_gaussian = exp(-pow(rad, kp['sharpness']))

            elif kp['type'] == 'above':
                # Place the cylinders above the rotor aligned with the Z-direction
                ax = (xs[2])/(0.5*kp['length'])      
                axial_gaussian = exp(-pow(ax, kp['sharpness']))
                rad = (xs[0]**2 + xs[1]**2)/kp['radius']**2
                radial_gaussian = exp(-pow(rad, kp['sharpness']))

            kernel_exp += axial_gaussian*radial_gaussian

            kernel_exp_list.append(axial_gaussian*radial_gaussian)

        return kernel_exp, kernel_exp_list

    # ================================================================

    # Modify some of the properties of the kernel
    kwargs['radius'] *= solver.problem.farm.turbines[0].RD
    kwargs['length'] *= solver.problem.farm.turbines[0].RD

    kernel_exp, kernel_exp_list = build_cylindrical_kernels(solver, kwargs)

    # Normalize this kernel function
    # vol = assemble(kernel_exp*dx)/solver.problem.farm.numturbs
    vol = assemble(kernel_exp*dx)
    kernel_exp = kernel_exp/vol

    save_debugging_files = False

    if save_debugging_files:
        folder_string = solver.params.folder+"data/"

        fp = File('%skernel_test.pvd' % (folder_string))
        kernel = project(kernel_exp, solver.problem.fs.Q, solver_type='cg')
        kernel.rename('kernel', 'kernel')
        fp << kernel

    # J = -assemble(sqrt(inner(solver.problem.u_k, solver.problem.u_k))*kernel_exp*dx)
    J = assemble(solver.problem.u_k[0]*kernel_exp*dx)


    calculate_individual_blockage = True

    print('Objective Value: ', float(J))

    if calculate_individual_blockage:

        blockage_array = []

        print('Calculating Individual Turbine Blockage')

        for k in range(solver.problem.farm.numturbs):
            kernel_exp_list[k] = kernel_exp_list[k]/vol

            # J_ind = -assemble(sqrt(inner(solver.problem.u_k, solver.problem.u_k))*kernel_exp_list[k]*dx)
            J_ind = assemble(solver.problem.u_k[0]*kernel_exp_list[k]/vol/solver.problem.farm.numturbs*dx)

            mx =  solver.problem.farm.turbines[k].mx
            my =  solver.problem.farm.turbines[k].my
            mz =  solver.problem.farm.turbines[k].mz
            yaw = solver.problem.farm.turbines[k].myaw

            blockage_array.append([k, mx, my, mz, yaw, J_ind])

        blockage_array = np.array(blockage_array)

        folder_string = solver.params.folder+"data/"
        np.savetxt('%sindividual_blockage.csv' % (folder_string),
            blockage_array,
            fmt='%.6e',
            header='Turbine ID (#), X-Location (m), Y-Location (m), Z-Location (m), Yaw (Rad), Blockage (m/s)',
            delimiter=',')

    return J
