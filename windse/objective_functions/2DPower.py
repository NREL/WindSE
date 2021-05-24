#######################################################################
####################### Preamble (do not edit) ########################
#######################################################################

import __main__
import os

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

#######################################################################
#######################################################################

### Additional import statements ###
import numpy as np
import math

### Declare Unique name
name = "2d_power"

### Define objective function
def objective(solver,inflow_angle = 0.0,first_call=False):
        x=SpatialCoordinate(solver.problem.dom.mesh)
        J=0.
        J_list=np.zeros(solver.problem.farm.numturbs+2)
        J_list[0]=solver.simTime
        for i in range(solver.problem.farm.numturbs):

            mx = solver.problem.farm.mx[i]
            my = solver.problem.farm.my[i]
            mz = solver.problem.farm.mz[i]
            x0 = [mx,my,mz]
            W = solver.problem.farm.W[i]*1.0
            R = solver.problem.farm.RD[i]/2.0 
            ma = solver.problem.farm.ma[i]
            yaw = solver.problem.farm.myaw[i]+delta_yaw
            u = solver.u_next
            A = pi*R**2.0
            C_tprime = 4*ma/(1-ma)
            C_pprime = 0.45/(1-ma)**3
            
            ### Rotate and Shift the Turbine ###
            xs = solver.problem.farm.YawTurbine(x,x0,yaw)
            u_d = u[0]*cos(yaw) + u[1]*sin(yaw)

            ### Create the function that represents the Thickness of the turbine ###
            T = exp(-pow((xs[0]/W),6.0))

            # WTGbase = Expression(("cos(yaw)","sin(yaw)"),yaw=yaw,degree=1)
            WTGbase = as_vector((cos(yaw),sin(yaw)))

            ### Create the function that represents the Disk of the turbine
            D = exp(-pow((pow((xs[1]/R),2)),6.0))

            volNormalization = assemble(T*D*dx)

            ### Create the function that represents the force ###
            if solver.problem.farm.force == "constant":
                F = 0.5*solver.problem.farm.RD[i]*C_tprime    
            elif solver.problem.farm.force == "sine":
                r = sqrt(xs[1]**2.0+xs[2]**2)
                F = 0.5*solver.problem.farm.RD[i]*C_tprime*(r/R*sin(pi*r/R)+0.5)/(.81831)

            J += (assemble(((0.5*A*C_pprime)**(1/3))*F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
            if solver.save_power or solver.save_objective:
                J_list[i+1] = 0.5*A*C_pprime*(assemble(F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
        
        if solver.save_power or solver.save_objective:
            J_list[-1]=float(J)

            folder_string = solver.params.folder+"data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)

            if first_call:
                f = open(folder_string+"2d_power_data.txt",'w')
                header = srt("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
                f.write(header)
            else:
                f = open(folder_string+"2d_power_data.txt",'a')

            np.savetxt(f,[J_list])
            f.close()

        return J