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
name = "KE_entrainment"

### Define objective function
def objective(solver,inflow_angle = 0.0,first_call=False):
    solver.fprint("Using Kinetic Energy Entrainment Functional")
    turb_id = solver.opt_turb_id[0]
    HH = solver.problem.farm.HH[turb_id]
    R = solver.problem.farm.RD[turb_id]/2.0

    # mark cells in area of interest
    if first_call:
        solver.KE_objective_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim())
        # solver.KE_objective_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim()-1)

        solver.KE_objective_markers.set_all(0)
        x0 = min(solver.problem.farm.x)
        x1 = solver.problem.dom.x_range[1]
        y0 = min(solver.problem.farm.y)-3.0*max(solver.problem.farm.RD)/2.0
        y1 = max(solver.problem.farm.y)+3.0*max(solver.problem.farm.RD)/2.0
        # if solver.ke_location =="tip":
        #     z0 = min(solver.problem.farm.z)
        #     z1 = max(solver.problem.farm.z)+max(solver.problem.farm.RD)/2.0
        # elif solver.ke_location == "hub":
        #     z0 = min(solver.problem.farm.z)-solver.problem.dom.mesh.hmin()*1.0
        #     z1 = max(solver.problem.farm.z)+solver.problem.dom.mesh.hmin()*1.0
        # elif solver.ke_location == "rotor":
        z0_in = min(solver.problem.farm.z) + R - solver.problem.dom.mesh.hmin()*1.0
        z1_in = max(solver.problem.farm.z) + R + solver.problem.dom.mesh.hmin()*1.0
        z0_out = min(solver.problem.farm.z) - R - solver.problem.dom.mesh.hmin()*1.0
        z1_out = max(solver.problem.farm.z) - R + solver.problem.dom.mesh.hmin()*1.0
 
        fluxIn =  CompiledSubDomain("x[0]>x0 && x[0]<x1 && x[1]>y0 && x[1]<y1  && x[2]>z0 && x[2]<z1",x0=x0,x1=x1,y0=y0,y1=y1,z0=z0_in,z1=z1_in)
        fluxOut = CompiledSubDomain("x[0]>x0 && x[0]<x1 && x[1]>y0 && x[1]<y1  && x[2]>z0 && x[2]<z1",x0=x0,x1=x1,y0=y0,y1=y1,z0=z0_out,z1=z1_out)
        # AOI  = CompiledSubDomain("x[0]>x0 && x[0]<x1 && x[1]>y0 && x[1]<y1  && x[2]>z0 && x[2]<z1",x0=x0,x1=x1,y0=y0,y1=y1,z0=z0,z1=z1)
        # AOI.mark(solver.KE_objective_markers,1)
        fluxIn.mark(solver.KE_objective_markers,1)
        fluxOut.mark(solver.KE_objective_markers,2)
        File(solver.params.folder+"mesh/KE_objective_markers.pvd")<<solver.KE_objective_markers

        if solver.save_objective:
            folder_string = solver.params.folder+"data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)
            f = open(folder_string+"KE_entrainment_data.txt",'w')
            header = str("Time    "+"Entrainment\n")
            f.write(header)
            f.close()

    dx_KE = Measure('dx', subdomain_data=solver.KE_objective_markers)
    # ds_KE = Measure('ds', subdomain_data=solver.KE_objective_markers)


    J = assemble(-(solver.problem.vertKE*dx_KE(1) - solver.problem.vertKE*dx_KE(2)))

    if solver.save_objective:
        folder_string = solver.params.folder+"data/"
        f = open(folder_string+"KE_entrainment_data.txt",'a')
        out_data = [solver.simTime]
        out_data.extend([float(J)])
        np.savetxt(f,[out_data])
        f.close()

    return J