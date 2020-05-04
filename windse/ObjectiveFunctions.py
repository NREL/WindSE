import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np
    import math


    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *



def CalculatePowerFunctional(solver,inflow_angle = 0.0):
    J = assemble(dot(solver.problem.tf,solver.problem.u_k)*dx)

    if solver.save_power:
        J_list=np.zeros(solver.problem.farm.numturbs+1)
        if solver.problem.farm.actuator_disks_list is not None:
            for i in range(solver.problem.farm.numturbs):
                yaw = solver.problem.farm.myaw[i]+inflow_angle
                tf1 = solver.problem.farm.actuator_disks_list[i] * cos(yaw)**2
                tf2 = solver.problem.farm.actuator_disks_list[i] * sin(yaw)**2
                tf3 = solver.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
                tf = tf1*solver.u_k[0]**2+tf2*solver.u_k[1]**2+tf3*solver.u_k[0]*solver.u_k[1]
                J_list[i] = assemble(dot(tf,solver.u_k)*dx,**solver.extra_kwarg)
        J_list[-1]=sum(J_list)

        folder_string = solver.params.folder+"/data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if solver.J_saved:
            f = open(folder_string+"power_data.txt",'ab')
        else:
            f = open(folder_string+"power_data.txt",'wb')
            solver.J_saved = True

        np.savetxt(f,[J_list])
        f.close()

    return J

def Calculate2DPowerFunctional(solver,inflow_angle = 0.0):
        x=SpatialCoordinate(solver.problem.dom.mesh)
        J=0.
        J_list=np.zeros(solver.problem.farm.numturbs+1)
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
            if solver.save_power:
                J_list[i] = 0.5*A*C_pprime*(assemble(F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
        
        if solver.save_power:
            J_list[-1]=np.sum(J_list[:-1])

            folder_string = solver.params.folder+"/data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)

            if solver.J_saved:
                f = open(folder_string+"power_data.txt",'ab')
            else:
                f = open(folder_string+"power_data.txt",'wb')
                solver.J_saved = True

            np.savetxt(f,[J_list])
            f.close()

        return J

def CalculateWakeCenter(solver,inflow_angle = 0.0):

    ### Get the maximum number of Roter Diameters down stream ###
    nRD = math.ceil((solver.problem.dom.x_range[1]-solver.problem.farm.x[0])/(solver.problem.farm.RD[0]))
    record_RD = min(nRD,solver.wake_RD)

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
        File(solver.params.folder+"/mesh/outflow_markers.pvd") << solver.outflow_markers

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
        if i+1 == record_RD:
            annotate = True

        ### Switch measure depending on location of RD ###
        if abs(x0[i] - solver.problem.dom.x_range[1]) <= 1e-2:
            print("External for " +repr(x0[i]))
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
    f.close()

    return J

def CalculateWakeCapture(solver,inflow_angle = 0.0):

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
            RDwindow  = CompiledSubDomain("near(x[0], x0, tol) && x[1]>=ly && x[1]<=uy  && x[2]>=lz && x[2]<=uz",x0 = x, ly=ly, uy=uy, lz=lz, uz=uz, tol = 1e-2)
            RDfull    = CompiledSubDomain("near(x[0], x0, tol)", x0 = x, tol = 1e-2)
            RDfull.mark(solver.outflow_markers,i+1)
            RDwindow.mark(solver.outflow_markers,2*(i+1))

        ### Save the markers for debugging ###
        File(solver.params.folder+"/mesh/outflow_markers.pvd") << solver.outflow_markers

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
    u_dif_mag = sqrt((u_ref[0]-u_ref[0])**2.0+(u_ref[1]-u_ref[1])**2.0+(u_ref[2]-u_ref[2])**2.0)
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
            # print("External for " +repr(x0[i]))
            M_full    = assemble(u_dif_mag*ds_external((i+1,2*(i+1))), annotate = annotate)
            M_window  = assemble(u_dif_mag*ds_external(2*(i+1)), annotate = annotate)
        else:
            M_full    = assemble(avg(u_dif_mag)*ds_internal((i+1,2*(i+1))), annotate = annotate)
            M_window  = assemble(avg(u_dif_mag)*ds_internal(2*(i+1)), annotate = annotate)
        M = M_window/M_full

        ### Collect Data ###
        out_data.extend([M])

        ### Return Objective Function ###
        if annotate:
            print("RD"+repr(i+1)+" Capture %: "+repr(M*100))
            J = M

    ### Save Data ###
    np.savetxt(f,[out_data])
    f.close()

    return J

objectives_dict = {"power":    CalculatePowerFunctional,
                   "2d-power":    Calculate2DPowerFunctional,
                   "wake_deflection": CalculateWakeCenter,
                   "wake_capture": CalculateWakeCapture
                  }