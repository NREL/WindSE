import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *



def CalculatePowerFunctional(solver,inflow_angle = 0.0):
    J = assemble(dot(solver.problem.tf,solver.problem.u_k)*dx)
    return J

def CalculateWakeCenter(solver,inflow_angle = 0.0):

    ### Get the maximum number of Roter Diameters down stream ###
    nRD = int(round((solver.problem.dom.x_range[1]-solver.problem.farm.x[0])/(solver.problem.farm.RD[0])))

    ### Get the mesh nodes that are closest to each Rotor Diameter ###
    xunique = np.unique(solver.problem.dom.mesh.coordinates()[:,0])
    x0 = []
    for i in range(nRD):
        xtarget = solver.problem.farm.x[0]+(i+1)*solver.problem.farm.RD[0]
        x0.append(xunique[np.argmin(np.abs(xunique-xtarget))])

    ### If we haven't made the markers yet, do it ###
    if not hasattr(solver,"outflow_markers"):

        ### Get some parameters ###
        y_factor = 2.0
        z_factor = 1.2
        HH = solver.problem.farm.HH[0]
        y0 = solver.problem.farm.y[0]
        R = solver.problem.farm.RD[0]/2.0
        ly=y0-y_factor*R
        uy=y0+y_factor*R
        lz=HH-z_factor*R
        uz=HH+z_factor*R

        ### Create the Facet Function ###
        solver.outflow_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim() - 1)
        solver.outflow_markers.set_all(0)

        ### Mark the mesh ###
        for i, x in enumerate(x0):
            RDregion  = CompiledSubDomain("near(x[0], x0, tol) && x[1]>=ly && x[1]<=uy  && x[2]>=lz && x[2]<=uz",x0 = x, ly=ly, uy=uy, lz=lz, uz=uz, tol = 1e-2)
            RDregion.mark(solver.outflow_markers,i+1)

        ### Save the markers for debugging ###
        # File("test.pvd") << solver.outflow_markers

        ### Create output folder ###
        folder_string = solver.params.folder+"/data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        ### Open the save file ### 
        f = open(folder_string+"wake_data.txt",'wb')

    else:

        ### Open to save file (append mode) ###
        folder_string = solver.params.folder+"/data/"
        f = open(folder_string+"wake_data.txt",'ab')

    ### Create the measures ###
    ds_internal = Measure('dS', subdomain_data=solver.outflow_markers)
    ds_external = Measure('ds', subdomain_data=solver.outflow_markers)
    x = SpatialCoordinate(solver.problem.dom.mesh)

    ### Get the 'Mass' Function ###
    u_ref = solver.problem.bd.bc_velocity
    u     = solver.problem.u_k
    u_dif_mag = sqrt((u[0]-u_ref[0])**2.0+(u[1]-u_ref[1])**2.0+(u[2]-u_ref[2])**2.0)
    # u_dif_mag = sqrt((u[0]-u_ref[0])**2.0)

    ### Calculate the Centroids for each RD ###
    out_data = [solver.simTime]
    for i in range(nRD):

        ### Check if this is the RD we want sensitivities ###
        annotate = False
        if i+1 == solver.wake_RD:
            annotate = True

        ### Switch measure depending on location of RD ###
        if abs(x0[i] - solver.problem.dom.x_range[1]) <= 1e-2:
            M = assemble(u_dif_mag*ds_external(i+1), annotate = annotate)
            Mx = assemble(x[0]*u_dif_mag*ds_external(i+1), annotate = False)
            My = assemble(x[1]*u_dif_mag*ds_external(i+1), annotate = annotate)
            Mz = assemble(x[2]*u_dif_mag*ds_external(i+1), annotate = False)
        else:
            M = assemble(avg(u_dif_mag)*ds_internal(i+1), annotate = annotate)
            Mx = assemble(avg(x[0]*u_dif_mag)*ds_internal(i+1), annotate = False)
            My = assemble(avg(x[1]*u_dif_mag)*ds_internal(i+1), annotate = annotate)
            Mz = assemble(avg(x[2]*u_dif_mag)*ds_internal(i+1), annotate = False)

        ### Collect Data ###
        out_data.extend([Mx/M,My/M,Mz/M])

        ### Return Objective Function ###
        if annotate:
            print("RD"+repr(i+1)+" Centroid: "+repr((Mx/M,My/M,Mz/M)))
            J = My/M

    ### Save Data ###
    np.savetxt(f,[out_data])

    return J

objectives_dict = {"power":    CalculatePowerFunctional,
                   "wake_deflection": CalculateWakeCenter
                  }