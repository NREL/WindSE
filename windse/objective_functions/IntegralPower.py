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
import math
import os

### Declare Unique name
name = "power"

### Set default keyword argument values ###
keyword_defaults = {}

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, annotate=True, **kwargs):
    '''
    The "power" objective function calculates the power using actuator disks, 
    by computing the integral of turbine force dotted with velocity. 
    '''

    # if not annotate:
    #     stop_annotating()

    J = assemble(dot(-solver.problem.tf,solver.problem.u_k)*dx)
    # J = assemble(dot(-solver.problem.tf,Constant((1.0,1.0,1.0)))*dx)

    # coords = solver.problem.fs.V.tabulate_dof_coordinates()
    # coords = coords[0::3, :]

    # ymin = solver.problem.dom.y_range[0]
    # ymax = solver.problem.dom.y_range[1]
    # y_coords = coords[:, 1]

    # zero_field = Function(solver.problem.fs.V)
    # temp_vector = np.copy(zero_field.vector()[:])
    # temp_vector[0::3] = -(y_coords - ymin)*(y_coords - ymax)

    # temp_vector = temp_vector/np.amax(temp_vector)*solver.problem.params['boundary_conditions']['HH_vel']

    # zero_field.vector()[:] = temp_vector[:]

    # saved_zero_field = File('zero_field.pvd')
    # saved_zero_field << zero_field

    # if solver.params['solver']['type'] == 'unsteady':
    #     fac = -1.0
    # else:
    #     fac = 1.0

    # J = -assemble(dot(fac*solver.problem.tf,zero_field)*dx)

    if solver.iter_val is None:
        iter_val = -1
    else:
        iter_val = solver.iter_val

    if solver.save_power or solver.save_objective:
        J_list=np.zeros(solver.problem.farm.numturbs+3)
        J_list[0]=iter_val
        J_list[1]=solver.simTime
        if getattr(solver.problem.farm,"actuator_disks_list",None) is not None:
            for i in range(solver.problem.farm.numturbs):
                yaw = solver.problem.farm.myaw[i]+inflow_angle
                tf1 = solver.problem.farm.actuator_disks_list[i] * cos(yaw)**2
                tf2 = solver.problem.farm.actuator_disks_list[i] * sin(yaw)**2
                tf3 = solver.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
                tf = tf1*solver.problem.u_k[0]**2+tf2*solver.problem.u_k[1]**2+tf3*solver.problem.u_k[0]*solver.problem.u_k[1]
                J_list[i+2] = assemble(dot(-tf,solver.problem.u_k)*dx,**solver.extra_kwarg)
            else:
                pass
                # print("WARNING: missing individual turbine actuator disk, only able to report full farm power")

        J_list[-1]=float(J)

        folder_string = solver.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if first_call:
            f = open(folder_string+"power_data.txt",'w')
            header = str("Iter_val    "+"Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
            f.write(header)
        else:
            f = open(folder_string+"power_data.txt",'a')

        np.savetxt(f,[J_list])
        f.close()

    # start_annotating()
    return J
