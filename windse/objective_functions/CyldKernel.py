### These must be imported ###
from dolfin import *
from dolfin_adjoint import *

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
    def rotate_and_shift_points(x, x0, yaw, kp):

        xs =  cos(yaw)*(x[0]-x0[0]) + sin(yaw)*(x[1]-x0[1])
        ys = -sin(yaw)*(x[0]-x0[0]) + cos(yaw)*(x[1]-x0[1])

        if x.geometric_dimension() == 3:
            zs = x[2] - x0[2]
        else:
            zs = 0.0

        if kp['type'] == 'upstream':
            xs = xs + 0.5*kp['length']
        elif kp['type'] == 'above':
            zs = zs - 0.5*kp['length']

        return [xs, ys, zs]

    def build_cylindrical_kernels(solver, kp):
        x = SpatialCoordinate(solver.problem.dom.mesh)

        kernel_exp = 0

        for i in range(solver.problem.farm.numturbs):

            # Get a convenience copy of turbine i's location
            mx = solver.problem.farm.mx[i]
            my = solver.problem.farm.my[i]
            mz = solver.problem.farm.mz[i]
            x0 = [mx, my, mz]

            # Get a convenience copy of turbine i's yaw
            yaw = solver.problem.farm.myaw[i]

            xs = rotate_and_shift_points(x, x0, yaw, kp)

            if kp['type'] == 'upstream':
                # Place the cylinders upstream from the rotor aligned with the hub axis
                ax = xs[0]/(0.5*kp['length'])
                axial_gaussian = exp(-pow(ax, kp['sharpness']))
                rad = (xs[1]**2 + xs[2]**2)/kp['radius']**2
                radial_gaussian = exp(-pow(rad, kp['sharpness']))

            elif kp['type'] == 'above':
                # Place the cylinders above the rotor aligned with the Z-direction
                ax = xs[2]/(0.5*kp['length'])        
                axial_gaussian = exp(-pow(ax, kp['sharpness']))
                rad = (xs[0]**2 + xs[1]**2)/kp['radius']**2
                radial_gaussian = exp(-pow(rad, kp['sharpness']))

            kernel_exp += axial_gaussian*radial_gaussian

        return kernel_exp

    # ================================================================

    # Modify some of the properties of the kernel
    kwargs['radius'] *= solver.problem.farm.RD[0]
    kwargs['length'] *= solver.problem.farm.RD[0]

    kernel_exp = build_cylindrical_kernels(solver, kwargs)

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

    return J
