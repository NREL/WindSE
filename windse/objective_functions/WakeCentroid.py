################## HEADER DO NOT EDIT ##################
import os
import __main__

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    from dolfin_adjoint import *
########################################################

### Additional import statements ###
import numpy as np
import math

### Declare Unique name
name = "wake_center"

### Set default keyword argument values ###
keyword_defaults = {
    "wake_RD":     5,
    "wake_length": 5,
    "wake_radius": 0.5
    }

### Define objective function
def objective(solver, inflow_angle = 0.0, first_call=False, **kwargs):
    '''
    The "wake_center" objective function computes the wake n rotor diameters
    downstream from a single turbine. This calculation is perform by centering
    a cylinder oriented along the streamwise direction downstream and calculating
    the center of mass of the velocity deficit, or centroid. 

    Keyword arguments:
        wake_RD:     Number of rotor diameters downstream where the centroid will be computed
        wake_length: The streamwise length for the area of integration (not used)
        wake_radius: The radius of the cylinder (not used)
    '''

    ### Extract keyword arguments
    wake_RD = kwargs.pop("wake_RD")
    wake_length = kwargs.pop("wake_length")
    wake_radius = kwargs.pop("wake_radius")


    turb_id = solver.opt_turb_id[0]

    ### Get the maximum number of Roter Diameters down stream ###
    nRD = math.ceil((solver.problem.dom.x_range[1]-solver.problem.farm.x[turb_id])/(solver.problem.farm.RD[turb_id]))
    record_RD = min(nRD,int(wake_RD))

    ### Get the mesh nodes that are closest to each Rotor Diameter ###
    xunique = np.unique(solver.problem.dom.mesh.coordinates()[:,0])
    x0 = []
    for i in range(nRD):
        xtarget = solver.problem.farm.x[turb_id]+(i+1)*solver.problem.farm.RD[turb_id]
        x0.append(xunique[np.argmin(np.abs(xunique-xtarget))])

    ### If we haven't made the markers yet, do it ###
    if first_call:

        ### Get some parameters ###
        y_factor = 2.0
        z_factor = 1.2
        HH = solver.problem.farm.HH[turb_id]
        y0 = solver.problem.farm.y[turb_id]
        R = solver.problem.farm.RD[turb_id]/2.0
        ly=y0-y_factor*R
        uy=y0+y_factor*R
        lz=HH-z_factor*R
        uz=HH+z_factor*R

        ### Create the Facet Function ###
        solver.WC_objective_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim() - 1)
        solver.WC_objective_markers.set_all(0)

        ### Mark the mesh ###
        for i, x in enumerate(x0):
            RDregion  = CompiledSubDomain("near(x[0], x0, tol) && x[1]>=ly && x[1]<=uy  && x[2]>=lz && x[2]<=uz",x0 = x, ly=ly, uy=uy, lz=lz, uz=uz, tol = 1e-2)
            RDregion.mark(solver.WC_objective_markers,i+1)

        if solver.save_objective:
            ### Save the markers for debugging ###
            File(solver.params.folder+"mesh/WC_objective_markers.pvd") << solver.WC_objective_markers

            ### Create output folder ###
            folder_string = solver.params.folder+"data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)

            ### Open the save file ### 
            f = open(folder_string+"wake_center_data.txt",'w')
            header = str("Time    "+"%dRD_cx    %dRD_cy    %dRD_cz    "*nRD % tuple(np.repeat(range(1,nRD+1),3))+"\n")
            f.write(header)

    else:
        if solver.save_objective:
            ### Open to save file (append mode) ###
            folder_string = solver.params.folder+"data/"
            f = open(folder_string+"wake_center_data.txt",'a')

    ### Create the measures ###
    ds_internal = Measure('dS', subdomain_data=solver.WC_objective_markers)
    ds_external = Measure('ds', subdomain_data=solver.WC_objective_markers)
    x = SpatialCoordinate(solver.problem.dom.mesh)

    ### Get the 'Mass' Function ###
    u_ref = solver.problem.bd.bc_velocity
    u     = solver.problem.u_k
    u_dif_mag = sqrt((u[0]-u_ref[0])**2.0+(u[1]-u_ref[1])**2.0+(u[2]-u_ref[2])**2.0)
    # solver.fprint(type(u_dif_mag))
    # u_dif_mag = sqrt((u[0]-u_ref[0])**2.0)

    ### Calculate the Centroids for each RD ###
    out_data = [solver.simTime]
    for i in range(nRD):

        ### Check if this is the RD we want sensitivities ###
        anno_switch={}
        anno_false={}
        if solver.params.dolfin_adjoint:
            anno_false["annotate"] = False
            anno_switch["annotate"] = False
            if i+1 == record_RD:
                anno_switch["annotate"] = True

        ### Switch measure depending on location of RD ###
        if abs(x0[i] - solver.problem.dom.x_range[1]) <= 1e-2:
            # solver.fprint("External for " +repr(x0[i]))
            M = assemble(u_dif_mag*ds_external(i+1), **anno_switch)
            Mx = assemble(x[0]*u_dif_mag*ds_external(i+1), **anno_false)
            My = assemble(x[1]*u_dif_mag*ds_external(i+1), **anno_switch)
            Mz = assemble(x[2]*u_dif_mag*ds_external(i+1), **anno_false)
        else:
            M = assemble(avg(u_dif_mag)*ds_internal(i+1), **anno_switch)
            Mx = assemble(avg(x[0]*u_dif_mag)*ds_internal(i+1), **anno_false)
            My = assemble(avg(x[1]*u_dif_mag)*ds_internal(i+1), **anno_switch)
            Mz = assemble(avg(x[2]*u_dif_mag)*ds_internal(i+1), **anno_false)

        ### Collect Data ###
        out_data.extend([Mx/M,My/M,Mz/M])

        ### Return Objective Function ###
        if i+1 == record_RD:
            solver.fprint("RD"+repr(i+1)+" Centroid: "+repr((Mx/M,My/M,Mz/M)))
            # J = -pow(My/M*My/M,1.0/2.0) ## if need to be strictly positive
            J = My/M


    ### Save Data ###
    if solver.save_objective:
        np.savetxt(f,[out_data])
        f.close()

    return J