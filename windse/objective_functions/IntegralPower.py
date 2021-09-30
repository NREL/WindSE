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
    spherical_gaussian = Expression("exp(-pow(((pow((x[0] - x0),2) + pow((x[1] - y0),2) + pow((x[2] - z0),2))/r),6.0))",x0=0,y0=0,z0=100,r=500, degree=5)

    print("tf1: ",assemble(dot(solver.problem.tf,as_vector((1.0,0.0,0.0)))*dx))
    print("tf2: ",assemble(dot(solver.problem.tf,as_vector((0.0,1.0,0.0)))*dx))
    print("tf3: ",assemble(dot(solver.problem.tf,as_vector((0.0,0.0,1.0)))*dx))
    print("uk1: ",assemble(dot(spherical_gaussian*solver.problem.u_k,as_vector((1.0,0.0,0.0)))*dx))
    print("uk2: ",assemble(dot(spherical_gaussian*solver.problem.u_k,as_vector((0.0,1.0,0.0)))*dx))
    print("uk3: ",assemble(dot(spherical_gaussian*solver.problem.u_k,as_vector((0.0,0.0,1.0)))*dx))


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

    if solver.save_power or solver.save_objective:
        J_list=np.zeros(solver.problem.farm.numturbs+2)
        J_list[0]=solver.simTime
        if getattr(solver.problem.farm,"actuator_disks_list",None) is not None:
            for i in range(solver.problem.farm.numturbs):
                yaw = solver.problem.farm.myaw[i]+inflow_angle
                tf1 = solver.problem.farm.actuator_disks_list[i] * cos(yaw)**2
                tf2 = solver.problem.farm.actuator_disks_list[i] * sin(yaw)**2
                tf3 = solver.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
                tf = tf1*solver.problem.u_k[0]**2+tf2*solver.problem.u_k[1]**2+tf3*solver.problem.u_k[0]*solver.problem.u_k[1]
                J_list[i+1] = assemble(dot(-tf,solver.problem.u_k)*dx,**solver.extra_kwarg)
            else:
                print("WARNING: missing individual turbine actuator disk, only able to report full farm power")

        J_list[-1]=float(J)

        folder_string = solver.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if first_call:
            f = open(folder_string+"power_data.txt",'w')
            header = str("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
            f.write(header)
        else:
            f = open(folder_string+"power_data.txt",'a')

        np.savetxt(f,[J_list])
        f.close()

    # start_annotating()
    return J
