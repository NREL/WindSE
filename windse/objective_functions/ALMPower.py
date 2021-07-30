### These must be imported ###
from dolfin import *
from dolfin_adjoint import *

### Additional import statements ###
import numpy as np
import math
import os

### Declare Unique name
name = "alm_power"

### Set default keyword argument values ###
keyword_defaults = {
    "alm_power_type": "real" 
    }

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    The "alm_power" objective function computes the power using actuator lines 
    by dotting the turbine force in the rotor plane with the moment arm of the 
    turbine blade multiplied by angular velocity. Can be used for multiple
    turbines.

    Keyword arguments:
        alm_power_type: real or fake; 
                            real - dotting the turbine force with the moment arm of the turbine blade 
                            fake - simply multiply the turbine force by the velocity
    '''


    # J = assemble(solver.problem.tf_list[0]*dx)
    # J = assemble(solver.problem.tf_list[0][0]*solver.problem.tf_list[0][0]*dx)
    # J = assemble(-2.0*np.pi*solver.problem.rpm/60.0*abs(inner(solver.problem.tf_list[0],solver.problem.cyld_expr_list[0]))*dx)
    # J = assemble(-inner(solver.problem.tf_list[0],solver.problem.cyld_expr_list[0])**2.0*dx)

    # J = assemble(inner(solver.problem.u_k,solver.problem.u_k)*dx)
    # J = assemble(inner(solver.problem.tf,solver.problem.tf)*dx)
    # J = assemble(1e-6*2.0*np.pi*solver.problem.rpm/60.0*dot(solver.problem.tf,solver.problem.u_k)*dx)

    # cyld_expr = Expression(('sin(yaw)*(x[2]-zs)', '-cos(yaw)*(x[2]-zs)', '(x[1]-ys)*cos(yaw)-(x[0]-xs)*sin(yaw)'),
    #     degree=6,
    #     yaw=solver.problem.farm.yaw[0],
    #     xs=solver.problem.farm.x[0],
    #     ys=solver.problem.farm.y[0],
    #     zs=solver.problem.farm.z[0])

    # J = assemble(2.0*np.pi*solver.problem.rpm/60.0*dot(solver.problem.tf,cyld_expr)*dx)
    # J = assemble(inner(solver.problem.u_k,as_vector((1.0,0.0,0.0)))*dx)
    # J = assemble(inner(solver.problem.u_k,as_vector((1.0,0.0,0.0)))*dx)

    ### Extract keyword arguments
    alm_power_type = kwargs.pop("alm_power_type")

    J_list=np.zeros(solver.problem.farm.numturbs+2)
    J_list[0]=solver.simTime
    J = 0.0
    for i in range(solver.problem.farm.numturbs):
        if alm_power_type == "real":
            J_temp = assemble(1e-6*(2.0*np.pi*solver.problem.rpm/60.0)*inner(-solver.problem.tf_list[i], solver.problem.cyld_expr_list[i])*dx)
        elif alm_power_type == "fake":
            J_temp = assemble(1e-6*(2.0*np.pi*solver.problem.rpm/60.0)*dot(-solver.problem.tf_list[i],solver.problem.u_k)*dx)
        else:
            raise ValueError("Unknown ALM Power type: "+repr(alm_power_type))
        J_list[i+1] = J_temp
        J += J_temp
    J_list[-1] = float(J)


    if solver.save_power or solver.save_objective:

        folder_string = solver.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if first_call:
            f = open(folder_string+"alm_power_data.txt",'w')
            header = str("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
            f.write(header)
        else:
            f = open(folder_string+"alm_power_data.txt",'a')

        np.savetxt(f,[J_list])
        f.close()

    return J
