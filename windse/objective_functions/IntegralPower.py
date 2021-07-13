### These must be imported ###
from dolfin import *
from dolfin_adjoint import *

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

    J = assemble(dot(solver.problem.tf,solver.problem.u_k)*dx)

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
                J_list[i+1] = assemble(dot(tf,solver.problem.u_k)*dx,**solver.extra_kwarg)
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
