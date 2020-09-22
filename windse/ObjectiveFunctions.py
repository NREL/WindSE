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
    import numpy as np
    import math


    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *



def CalculatePowerFunctional(solver,inflow_angle = 0.0):
    J = -assemble(dot(solver.problem.tf,solver.problem.u_k)*dx)

    if solver.save_power or solver.save_objective:
        J_list=np.zeros(solver.problem.farm.numturbs+2)
        J_list[0]=solver.simTime
        if getattr(solver.problem.farm,"actuator_disks_list",None) is not None:
            for i in range(solver.problem.farm.numturbs):
                yaw = solver.problem.farm.myaw[i]+inflow_angle
                tf1 = solver.problem.farm.actuator_disks_list[i] * cos(yaw)**2
                tf2 = solver.problem.farm.actuator_disks_list[i] * sin(yaw)**2
                tf3 = solver.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
                tf = tf1*solver.u_k[0]**2+tf2*solver.u_k[1]**2+tf3*solver.u_k[0]*solver.u_k[1]
                J_list[i+1] = assemble(dot(tf,solver.u_k)*dx,**solver.extra_kwarg)
            else:
                print("WARNING: missing individual turbine actuator disk, only able to report full farm power")

        J_list[-1]=float(J)

        folder_string = solver.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if solver.pwr_saved:
            f = open(folder_string+"power_data.txt",'a')
        else:
            f = open(folder_string+"power_data.txt",'w')
            header = str("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
            f.write(header)
            solver.pwr_saved = True

        np.savetxt(f,[J_list])
        f.close()

    return J

def CalculateActuatorLinePowerFunctional(solver,inflow_angle = 0.0):
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

    J_list=np.zeros(solver.problem.farm.numturbs+2)
    J_list[0]=solver.simTime
    J = 0.0
    for i in range(solver.problem.farm.numturbs):
        if solver.alm_power_type == "real":
            J_temp = -assemble(1e-6*(2.0*np.pi*solver.problem.rpm/60.0)*inner(-solver.problem.tf_list[i], solver.problem.cyld_expr_list[i])*dx, annotate=True)
        elif solver.alm_power_type == "fake":
            J_temp = -assemble(1e-6*(2.0*np.pi*solver.problem.rpm/60.0)*dot(-solver.problem.tf_list[i],solver.problem.u_k)*dx)
        else:
            raise ValueError("Unknown ALM Power type: "+repr(solver.alm_power_type))
        J_list[i+1] = J_temp
        J += J_temp
    J_list[-1] = float(J)


    if solver.save_power or solver.save_objective:

        folder_string = solver.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if solver.pwr_saved:
            f = open(folder_string+"power_data.txt",'a')
        else:
            f = open(folder_string+"power_data.txt",'w')
            header = str("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
            f.write(header)
            solver.pwr_saved = True

        np.savetxt(f,[J_list])
        f.close()

    return J

def Calculate2DPowerFunctional(solver,inflow_angle = 0.0):
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

            if solver.pwr_saved:
                f = open(folder_string+"power_data.txt",'a')
            else:
                f = open(folder_string+"power_data.txt",'w')
                header = srt("Time    "+"Turbine_%d    "*solver.problem.farm.numturbs % tuple(range(solver.problem.farm.numturbs))+"Sum"+"\n")
                f.write(header)
                solver.pwr_saved = True

            np.savetxt(f,[J_list])
            f.close()

        return J

def CalculateKEEntrainment(solver,inflow_angle = 0.0):
    print("Using Kinetic Energy Entrainment Functional")

    # mark cells in area of interest
    if not hasattr(solver,"outflow_markers"):
        solver.objective_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim())
        solver.objective_markers.set_all(0)
        x0 = min(solver.problem.farm.x)
        x1 = solver.problem.dom.x_range[1]
        y0 = min(solver.problem.farm.y)-3.0*max(solver.problem.farm.RD)/2.0
        y1 = max(solver.problem.farm.y)+3.0*max(solver.problem.farm.RD)/2.0
        if solver.ke_location =="tip":
            z0 = min(solver.problem.farm.z)
            z1 = max(solver.problem.farm.z)+max(solver.problem.farm.RD)/2.0
        elif solver.ke_location == "hub":
            z0 = min(solver.problem.farm.z)-solver.problem.dom.mesh.hmin()*1.0
            z1 = max(solver.problem.farm.z)+solver.problem.dom.mesh.hmin()*1.0

        AOI  = CompiledSubDomain("x[0]>x0 && x[0]<x1 && x[1]>y0 && x[1]<y1  && x[2]>z0 && x[2]<z1",x0=x0,x1=x1,y0=y0,y1=y1,z0=z0,z1=z1)
        AOI.mark(solver.objective_markers,1)
        File(solver.params.folder+"test2.pvd")<<solver.objective_markers

        if solver.save_objective:
            folder_string = solver.params.folder+"data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)
            f = open(folder_string+"objective_data.txt",'w')
            header = str("Time    "+"Entrainment\n")
            f.write(header)
            f.close()

    dx_KE = Measure('dx', subdomain_data=solver.objective_markers)

    J = assemble(-1e-6*solver.problem.vertKE*dx_KE(1))

    if solver.save_objective:
        folder_string = solver.params.folder+"data/"
        f = open(folder_string+"objective_data.txt",'a')
        out_data = [solver.simTime]
        out_data.extend([float(J)])
        np.savetxt(f,[out_data])
        f.close()

    return J



def CalculateWakeCenter(solver,inflow_angle = 0.0):

    turb_id = solver.opt_turb_id[0]

    ### Get the maximum number of Roter Diameters down stream ###
    nRD = math.ceil((solver.problem.dom.x_range[1]-solver.problem.farm.x[turb_id])/(solver.problem.farm.RD[turb_id]))
    record_RD = min(nRD,solver.wake_RD)

    ### Get the mesh nodes that are closest to each Rotor Diameter ###
    xunique = np.unique(solver.problem.dom.mesh.coordinates()[:,0])
    x0 = []
    for i in range(nRD):
        xtarget = solver.problem.farm.x[turb_id]+(i+1)*solver.problem.farm.RD[turb_id]
        x0.append(xunique[np.argmin(np.abs(xunique-xtarget))])

    ### If we haven't made the markers yet, do it ###
    if not hasattr(solver,"outflow_markers"):

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
        solver.outflow_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim() - 1)
        solver.outflow_markers.set_all(0)

        ### Mark the mesh ###
        for i, x in enumerate(x0):
            RDregion  = CompiledSubDomain("near(x[0], x0, tol) && x[1]>=ly && x[1]<=uy  && x[2]>=lz && x[2]<=uz",x0 = x, ly=ly, uy=uy, lz=lz, uz=uz, tol = 1e-2)
            RDregion.mark(solver.outflow_markers,i+1)

        if solver.save_objective:
            ### Save the markers for debugging ###
            File(solver.params.folder+"mesh/outflow_markers.pvd") << solver.outflow_markers

            ### Create output folder ###
            folder_string = solver.params.folder+"data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)

            ### Open the save file ### 
            f = open(folder_string+"objective_data.txt",'w')
            header = str("Time    "+"%dRD_cx    %dRD_cy    %dRD_cz    "*nRD % tuple(np.repeat(range(1,nRD+1),3))+"\n")
            f.write(header)

    else:
        if solver.save_objective:
            ### Open to save file (append mode) ###
            folder_string = solver.params.folder+"data/"
            f = open(folder_string+"objective_data.txt",'a')

    ### Create the measures ###
    ds_internal = Measure('dS', subdomain_data=solver.outflow_markers)
    ds_external = Measure('ds', subdomain_data=solver.outflow_markers)
    x = SpatialCoordinate(solver.problem.dom.mesh)

    ### Get the 'Mass' Function ###
    u_ref = solver.problem.bd.bc_velocity
    u     = solver.problem.u_k
    u_dif_mag = sqrt((u[0]-u_ref[0])**2.0+(u[1]-u_ref[1])**2.0+(u[2]-u_ref[2])**2.0)
    print(type(u_dif_mag))
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
            # print("External for " +repr(x0[i]))
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
            # J = -pow(My/M*My/M,1.0/2.0) ## if need to be strictly positive
            J = -My/M


    ### Save Data ###
    if solver.save_objective:
        np.savetxt(f,[out_data])
        f.close()

    return J

def CalculateWakeDeficit(solver,inflow_angle = 0.0):
    turb_id = solver.opt_turb_id[0]

    ### If we haven't made the markers yet, do it ###
    if not hasattr(solver,"outflow_markers"):

        ### Get some parameters ###
        z0 = solver.problem.farm.HH[turb_id]
        x0 = solver.problem.farm.x[turb_id]
        y0 = solver.problem.farm.y[turb_id]
        RD = solver.problem.farm.RD[turb_id]
        R = solver.wake_radius*RD
        L = solver.wake_length*RD
        lx=x0
        ux=x0+L

        ### Create the Facet Function ###
        solver.outflow_markers = MeshFunction("size_t", solver.problem.dom.mesh, solver.problem.dom.mesh.topology().dim() - 1)
        solver.outflow_markers.set_all(0)

        ### Mark the mesh ###
        wake_region  = CompiledSubDomain("x[0]>=lx && x[0]<=ux && x[1]*x[1]+x[2]*x[2] <= R*R", lx=lx, ux=ux, R=R, tol = 1e-2)
        wake_region.mark(solver.outflow_markers,1)

        unit = Function(solver.problem.fs.Q)
        unit.vector()[:] = 1.0
        dx_internal = Measure('dx', subdomain_data=solver.outflow_markers)
        solver.outflow_markers_vol = assemble(unit*dx_internal(1))

        if solver.save_objective:
            ### Save the markers for debugging ###
            File(solver.params.folder+"mesh/outflow_markers.pvd") << solver.outflow_markers

            ### Create output folder ###
            folder_string = solver.params.folder+"data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)

            ### Open the save file ### 
            f = open(folder_string+"objective_data.txt",'w')

    else:
        if solver.save_objective:
            ### Open to save file (append mode) ###
            folder_string = solver.params.folder+"data/"
            f = open(folder_string+"objective_data.txt",'a')

    ### Create the measures ###
    dx_internal = Measure('dx', subdomain_data=solver.outflow_markers)

    ### Get the Deficit Function ###
    u_ref = solver.problem.bd.bc_velocity
    u     = solver.problem.u_k
    ux_ref,uy_ref,uz_ref = u_ref.split(True)
    ux,uy,uz = u.split(True)
    # u_dif_mag = sqrt((u[0]-u_ref[0])**2.0+(u[1]-u_ref[1])**2.0+(u[2]-u_ref[2])**2.0)
    u_dif_mag = (ux_ref-ux)

    ### Calculate the Centroids for each RD ###
    out_data = [solver.simTime]

    ### Calculate the percent change ###
    J_deficit = assemble(u_dif_mag*dx_internal(1))/solver.outflow_markers_vol
    J_ref     = assemble(ux_ref*dx_internal(1))/solver.outflow_markers_vol
    J         = assemble(u_dif_mag/ux_ref*dx_internal(1))/solver.outflow_markers_vol

    ### Collect Data ###
    out_data.extend([J_deficit,J_ref,J])

    ### Return Objective Function ###
    print("Wake Deficit %: "+repr(J*100))

    if solver.save_objective:
        ### Save Data ###
        np.savetxt(f,[out_data])
        f.close()

    return J

objectives_dict = {"power":    CalculatePowerFunctional,
                   "alm_power": CalculateActuatorLinePowerFunctional,
                   "2d_power":    Calculate2DPowerFunctional,
                   "wake_center": CalculateWakeCenter,
                   "wake_deficit": CalculateWakeDeficit,
                   "KE_entrainment":    CalculateKEEntrainment
                  }