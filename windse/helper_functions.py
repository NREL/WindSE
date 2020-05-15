import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from windse import windse_parameters
    if windse_parameters.dolfin_adjoint:
        from dolfin import dx, File
        from dolfin_adjoint import Constant, Function, assemble
    else:
        from dolfin import Constant, Function, dx, assemble, File

    import numpy as np
    import time
    from scipy.special import gamma

def BaseHeight(x,y,ground,dx=0,dy=0):
    return Constant(ground(float(x),float(y),dx=dx,dy=dy))

def Transform(x,x0,HH,yaw,ground,dfd=None):

    ### Get the dimension ###
    dim = x.shape[0]

    ### Rotate and Shift the kernel ###
    if dfd is None:
        xhat = np.subtract.outer(x[0],x0[0])
        yhat = np.subtract.outer(x[1],x0[1])
        xrot =   np.cos(yaw)*xhat + np.sin(yaw)*yhat
        yrot = - np.sin(yaw)*xhat + np.cos(yaw)*yhat
        if dim == 3:
            z0 = ground(x0[0],x0[1])+HH
            zrot =   np.subtract.outer(x[2],z0)
    elif dfd == "x":
        xrot = - np.cos(yaw) 
        yrot =   np.sin(yaw) 
        if dim == 3:
            zrot = - ground(x0[0],x0[1],dx=1)
    elif dfd == "y":
        xrot = - np.sin(yaw)
        yrot = - np.cos(yaw)
        if dim == 3:
            zrot = - ground(x0[0],x0[1],dy=1)
    elif dfd == "yaw":
        xhat = np.subtract.outer(x[0],x0[0])
        yhat = np.subtract.outer(x[1],x0[1])
        xrot = - np.sin(yaw)*xhat + np.cos(yaw)*yhat
        yrot = - np.cos(yaw)*xhat - np.sin(yaw)*yhat
        if dim == 3:
            zrot =   0.0
    else:
        xrot = 0.0
        yrot = 0.0
        zrot = 0.0

    ### Adjust for 2D simulations ###
    if dim != 3:
        zrot = 0.0

    return [xrot,yrot,zrot]

def CalculateDiskTurbineForces(x,wind_farm,fs,dfd=None,save_actuators=False,sparse_ids=None,sparse_RDs=1.5,tfs=None):
    
    ### Collect the relevant turbine data ###
    x0 = np.array([wind_farm.mx,wind_farm.my,wind_farm.mz],dtype=float)
    yaw = np.array(wind_farm.myaw,dtype=float)+wind_farm.inflow_angle
    a = np.array(wind_farm.ma,dtype=float)
    HH = wind_farm.HH
    W = wind_farm.W
    R = wind_farm.RD/2.0
    dim, N = x.shape

    ### Set up some dim dependent values ###
    S_norm = (2.0+pi)/(2.0*pi)
    T_norm = 2.0*gamma(7.0/6.0)
    if dim == 3:
        A = np.pi*R**2.0 
        D_norm = np.pi*gamma(4.0/3.0)
    else:
        A = 2*R 
        D_norm = 2.0*gamma(7.0/6.0)
    volNormalization = T_norm*D_norm*W*R**(dim-1)

    ### Calculate relevant dofs that will be nonzero ###
    if sparse_ids is None:
        # print("recalc sparse")
        ### Create the bounding box for each turbine ###
        bounding_limit = [sparse_RDs*R]*3
        bl = x0-bounding_limit
        bu = x0+bounding_limit

        if dim == 3:
            ### Check all coordinate for which if above lower limits ###
            gt = np.greater.outer(x[0],bl[0])*np.greater.outer(x[1],bl[1])*np.greater.outer(x[2],bl[2])

            ### Check all coordinate for which if below upper limits ###
            lt = np.less.outer(x[0],bu[0])*np.less.outer(x[1],bu[1])*np.less.outer(x[2],bu[2])
        else:
            ### Check all coordinate for which if above lower limits ###
            gt = np.greater.outer(x[0],bl[0])*np.greater.outer(x[1],bl[1])

            ### Check all coordinate for which if below upper limits ###
            lt = np.less.outer(x[0],bu[0])*np.less.outer(x[1],bu[1])

        ### Check which coordinate are in bounding boxes ###
        sparse_ids = np.where(np.any(np.logical_and(gt,lt),axis=1))[0]

    ### Select sparse x values ###
    x=x[:,sparse_ids]

    ### Define Radial Force Functions ###
    if wind_farm.force == "constant":
        def RForce(r): return 1.0
        def dRForce(r,d_r): return 0.0
    elif wind_farm.force == "sine":
        def RForce(r): return (r*np.sin(np.pi*r)+0.5)/S_norm
        def dRForce(r,d_r): return (r*np.cos(np.pi*r)*(np.pi*d_r) + d_r*np.sin(np.pi*r))/S_norm
    else:
        ValueError("Unknown force type: "+wind_farm.force)

    if dfd is None:
        xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground)
        r = np.sqrt(np.power(xrot[1],2.0)+np.power(xrot[2],2.0))/R
        disks = (
                # Force
                4.*0.5*A*a/(1.-a)*RForce(r) * 
                # Disk Kernel
                np.exp(-(np.power(r,6.0)+np.power(xrot[0]/W,6.0)))/volNormalization
                )

        ### Rotate the force for the yawed turbine ###
        actuators_x = disks*np.cos(yaw)
        actuators_y = disks*np.sin(yaw)

        ### Create the normal ###
        n1 = np.cos(yaw)**2
        n2 = np.sin(yaw)**2
        n3 = 2.0*np.cos(yaw)*np.sin(yaw)

        if tfs is None:
            ### Initialize the output ###
            tf1  = Function(fs.tf_V)
            tf2  = Function(fs.tf_V)
            tf3  = Function(fs.tf_V)
            tfs = [tf1,tf2,tf3]

        ### Fill the output
        n=[n1,n2,n3]
        temp = np.zeros((N*dim))
        for i in range(3):
            temp[dim*(sparse_ids)+0] = -np.sum(actuators_x*n[i],axis=1)
            temp[dim*(sparse_ids)+1] = -np.sum(actuators_y*n[i],axis=1)
            tfs[i].vector()[:] = temp

    elif dfd in ["x","y"]:
        xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground)
        d_xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground,dfd=dfd)
        r = np.sqrt(np.power(xrot[1],2.0)+np.power(xrot[2],2.0))/R
        d_r = 0.5/R**2.0*(2.0*xrot[1]*d_xrot[1]+2.0*xrot[2]*d_xrot[2])/r
        d_disks = (
                  # # Force
                  4.*0.5*A*a/(1.-a) *
                  (
                    RForce(r) *
                    # # Derivative of Disk Kernel
                    -(6.0*np.power(r,5.0)*d_r+6.0*np.power(xrot[0]/W,5.0)*d_xrot[0]/W) +
                    # # Derivative of Force
                    dRForce(r,d_r)
                    # Disk Kernal
                  ) *
                  np.exp(-(np.power(r,6.0)+np.power(xrot[0]/W,6.0)))/volNormalization
                  )

        d_actuators_x = d_disks*np.cos(yaw)
        d_actuators_y = d_disks*np.sin(yaw)

        ### Create the normal ###
        n1 = np.cos(yaw)**2
        n2 = np.sin(yaw)**2
        n3 = 2.0*np.cos(yaw)*np.sin(yaw)

        ### Initialize the output ###
        tf1  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))
        tf2  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))
        tf3  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))

        ### Fill the output
        tfs = [tf1,tf2,tf3]
        n=[n1,n2,n3]
        for i in range(3):
            tfs[i][0::dim] = -d_actuators_x*n[i]
            tfs[i][1::dim] = -d_actuators_y*n[i]

    elif dfd == "a":
        xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground)
        r = np.sqrt(np.power(xrot[1],2.0)+np.power(xrot[2],2.0))/R
        d_disks = (
                  # Derivative of Force
                  4*0.5*A/(a-1.)**2.0*RForce(r) *
                  # Disk Kernal
                  np.exp(-(np.power(r,6.0)+np.power(xrot[0]/W,6.0)))/volNormalization
                  )

        d_actuators_x = d_disks*np.cos(yaw)
        d_actuators_y = d_disks*np.sin(yaw)

        ### Create the normal ###
        n1 = np.cos(yaw)**2
        n2 = np.sin(yaw)**2
        n3 = 2.0*np.cos(yaw)*np.sin(yaw)

        ### Initialize the output ###
        tf1  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))
        tf2  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))
        tf3  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))

        ### Fill the output
        tfs = [tf1,tf2,tf3]
        n=[n1,n2,n3]
        for i in range(3):
            tfs[i][0::dim] = -d_actuators_x*n[i]
            tfs[i][1::dim] = -d_actuators_y*n[i]

    elif dfd == "yaw":
        xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground)
        d_xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground,dfd=dfd)
        r = np.sqrt(np.power(xrot[1],2.0)+np.power(xrot[2],2.0))/R
        d_r = 0.5/R**2.0*(2.0*xrot[1]*d_xrot[1]+2.0*xrot[2]*d_xrot[2])/r
        disks = (
                # Force
                4.*0.5*A*a/(1.-a)*RForce(r) * 
                # Disk Kernel
                np.exp(-(np.power(r,6.0)+np.power(xrot[0]/W,6.0)))/volNormalization
                )
        d_disks = (
                  # # Force
                  4.*0.5*A*a/(1.-a) *
                  (
                    RForce(r) *
                    # # Derivative of Disk Kernel
                    -(6.0*np.power(r,5.0)*d_r+6.0*np.power(xrot[0]/W,5.0)*d_xrot[0]/W) +
                    # # Derivative of Force
                    dRForce(r,d_r)
                    # Disk Kernal
                  ) *
                  np.exp(-(np.power(r,6.0)+np.power(xrot[0]/W,6.0)))/volNormalization
                  )

        actuators_x = disks*np.cos(yaw)
        actuators_y = disks*np.sin(yaw)
        d_actuators_x = d_disks*np.cos(yaw) - disks*np.sin(yaw)
        d_actuators_y = d_disks*np.sin(yaw) + disks*np.cos(yaw)

        ### Create the normal ###
        n1 = np.cos(yaw)**2
        n2 = np.sin(yaw)**2
        n3 = 2.0*np.cos(yaw)*np.sin(yaw)
        d_n1 = (-2)*np.cos(yaw)*np.sin(yaw)
        d_n2 = 2*np.sin(yaw)*np.cos(yaw)
        d_n3 = 2.0*(np.cos(2*yaw))

        ### Initialize the output ###
        tf1  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))
        tf2  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))
        tf3  = np.zeros((len(sparse_ids)*dim,wind_farm.numturbs))

        ### Fill the output
        tfs = [tf1,tf2,tf3]
        n=[n1,n2,n3]
        d_n=[d_n1,d_n2,d_n3]
        for i in range(3):
            tfs[i][0::dim] = -(d_actuators_x*n[i]+actuators_x*d_n[i])
            tfs[i][1::dim] = -(d_actuators_y*n[i]+actuators_y*d_n[i])

    else:
        raise ValueError("Cannot take the derivative with respect to: "+dfd)

    ### Output the actuator information if needed ###
    if save_actuators and dfd is None:
        actuator_array = np.zeros((N*dim,wind_farm.numturbs))
        actuator_array[dim*(sparse_ids)+0] = actuators_x
        actuator_array[dim*(sparse_ids)+1] = actuators_y
    else:
        actuator_array = None

    return [[tf1,tf2,tf3],sparse_ids,actuator_array]




#================================================================

def UpdateActuatorLineForce(problem, simTime, dfd):
    def rot_x(theta):
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(theta), -np.sin(theta)],
                       [0, np.sin(theta), np.cos(theta)]])

        return Rx

    def rot_y(theta):
        Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                       [0, 1, 0],
                       [-np.sin(theta), 0, np.cos(theta)]])
        
        return Ry

    def rot_z(theta):
        Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                       [np.sin(theta), np.cos(theta), 0],
                       [0, 0, 1]])
        
        return Rz

    def write_lift_and_drag(fn, simTime, theta, forceVec):
        fp = open('./output/%s/nodal_%s.csv' % (problem.params.name, fn), 'a')
        fp.write('%.5e, %.5e, ' % (simTime, theta))

        N = np.size(forceVec)

        for j in range(N):
            if j < N-1:
                fp.write('%.5e, ' % (forceVec[j]))
            else:
                fp.write('%.5e\n' % (forceVec[j]))
        fp.close()

    def save_derivative_file(filename, deriv_array):

        dolfin_function = Function(problem.fs.V)

        fp = File('output/timeSeries/%s.pvd' % (filename))

        for k in range(problem.num_blade_segments):
            vec = np.copy(deriv_array[:, k])

            # Zero-out all very small values (<1e-12)
            vec[np.abs(vec) < 1e-12] = 0.0

            # Set the values of the dolfin function being saved
            dolfin_function.vector()[:] = vec

            # Rename the function for ParaView naming consistency
            dolfin_function.rename(filename, filename)

            # Save the function            
            fp << (dolfin_function, k)


    #================================================================
    # Get Mesh Properties
    #================================================================

    ndim = problem.dom.dim

    # Get the coordinates of the vector function space
    coords = problem.fs.V.tabulate_dof_coordinates()
    coords = np.copy(coords[0::problem.dom.dim, :])


    # Resape a linear copy of the coordinates for every mesh point
    coordsLinear = np.copy(coords.reshape(-1, 1))

    #================================================================
    # Set Turbine and Fluid Properties
    #================================================================

    # Set the density
    rho = 1.0

    # Set the hub height
    hub_height = problem.farm.HH[0]

    # Set the yaw of the turbine
    yaw = problem.farm.yaw[0]

    # Set the number of blades in the turbine
    num_blades = 3

    # Blade length (turbine radius)
    L = problem.farm.radius[0]

    # Width of Gaussian
    # Note: this sets the gaussian width to roughly twice the minimum cell length scale
    eps = 2.0*problem.dom.mesh.hmin()/np.sqrt(3)

    #================================================================
    # Set Derived Constants
    #================================================================

    # Calculate the radial position of each actuator node
    rdim = np.linspace(0.0, L, problem.num_blade_segments)

    # Calculate width of an individual blade segment
    w = rdim[1] - rdim[0]

    # Calculate an array describing the x, y, z position of each actuator node
    # Note: The basic blade is oriented along the +y-axis
    blade_pos_base = np.vstack((np.zeros(problem.num_blade_segments),
                        rdim,
                        np.zeros(problem.num_blade_segments)))

    # Calculate the blade velocity
    tip_speed = np.pi*2.0*L*problem.rpm/60.0

    # Specify the velocity vector at each actuator node
    # Note: A blade with span oriented along the +y-axis moves in the +z direction
    blade_vel_base = np.vstack((np.zeros(problem.num_blade_segments),
                           np.zeros(problem.num_blade_segments),
                           np.linspace(0.0, tip_speed, problem.num_blade_segments)))

    # Set the spacing pf each blade
    theta_vec = np.linspace(0.0, 2.0*np.pi, num_blades, endpoint = False)

    # Create unit vectors aligned with blade geometry
    # blade_unit_vec_base[:, 0] = points along rotor shaft
    # blade_unit_vec_base[:, 1] = points along blade span axis
    # blade_unit_vec_base[:, 3] = points tangential to blade span axis (generates a torque about rotor shaft)
    blade_unit_vec_base = np.array([[1.0, 0.0, 0.0],
                               [0.0, 1.0, 0.0],
                               [0.0, 0.0, 1.0]])

    #================================================================
    # Begin Calculating Turbine Forces
    #================================================================

    # Lift and drag coefficient (could be an array and you interpolate the value based on R)
    # cl_dolf = Constant((np.linspace(1.5, 0.5, problem.num_blade_segments)))
    # cd_dolf = Constant((np.ones(problem.num_blade_segments)))
    # cl = cl_dolf.values()
    # cd = cd_dolf.values()

    # Read cl and cd from the values specified in problem manager
    cl = np.array(problem.mcl, dtype = float)
    cd = np.array(problem.mcd, dtype = float)

    # Read the chord length from the values specified in the problem manager
    # c = L/20.0
    c = np.array(problem.mchord, dtype = float)

    # Initialize the sum of the torque on the rotor shaft to zero
    problem.rotor_torque = 0.0

    # Initialze arrays depending on what this function will be returning
    if dfd is None:
        tf_vec = np.zeros(np.size(coords))
        lift_force = np.zeros((np.shape(coords)[0], ndim))
        drag_force = np.zeros((np.shape(coords)[0], ndim))

    elif dfd == 'c_lift':
        cl = np.ones(problem.num_blade_segments)
        dfd_c_lift = np.zeros((np.size(coords), problem.num_blade_segments))

    elif dfd == 'c_drag':
        cd = np.ones(problem.num_blade_segments)
        dfd_c_drag = np.zeros((np.size(coords), problem.num_blade_segments))

    elif dfd == 'chord':
        c = np.ones(problem.num_blade_segments)
        dfd_chord = np.zeros((np.size(coords), problem.num_blade_segments))

    save_debugging_files = True

    if save_debugging_files:
        blade_pos_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        lift_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        drag_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        torque_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))


    # Calculate the blade position based on current simTime and turbine RPM
    period = 60.0/problem.rpm
    theta_offset = simTime/period*2.0*np.pi


    # Treat each blade separately
    for blade_ct, theta_0 in enumerate(theta_vec):
        theta = theta_0 + theta_offset

        # Generate a rotation matrix for this turbine blade
        Rx = rot_x(theta)
        Rz = rot_z(yaw)

        # Rotate the entire [x; y; z] matrix using this matrix, then shift to the hub height
        blade_pos = np.dot(Rz, np.dot(Rx, blade_pos_base))
        blade_pos[2, :] += hub_height

        # Initialize space to hold the fluid velocity at each actuator node
        u_fluid = np.zeros((3, problem.num_blade_segments))

        # Set if using local velocity around inidividual nodes
        using_local_velocity = False

        if using_local_velocity:
            # Generate the fluid velocity from the actual node locations in the flow
            for k in range(problem.num_blade_segments):
                u_fluid[:, k] = problem.u_k1(blade_pos[0, k],
                                             blade_pos[1, k],
                                             blade_pos[2, k])
                                
        else:
            # Generate the fluid velocity using the inlet velocity (x-pos = x_range[0])        
            for k in range(problem.num_blade_segments):
                u_fluid[:, k] = problem.u_k1(problem.dom.x_range[0],
                                             blade_pos[1, k],
                                             blade_pos[2, k])
                
                # Force inlet velocity (1, 0, 0) for debugging
                # u_fluid[:, k] = np.array([1, 0, 0])

        
        # Rotate the blade velocity in the global x, y, z, coordinate system
        # Note: blade_vel_base is negative since we seek the velocity of the fluid relative to a stationary blade
        # and blade_vel_base is defined based on the movement of the blade
        blade_vel = np.dot(Rz, np.dot(Rx, -blade_vel_base))
                        
        # Form the total relative velocity vector (including velocity from rotating blade)
        u_rel = u_fluid + blade_vel

        # Force fluid-only relative velocity for debugging
        # u_rel = u_fluid

        # Create unit vectors in the direction of u_rel
        u_rel_mag = np.linalg.norm(u_rel, axis=0)
        u_rel_mag[u_rel_mag < 1e-6] = 1e-6
        u_unit_vec = u_rel/u_rel_mag
        
        # Calculate the lift and drag forces using the relative velocity magnitude
        lift = (0.5*cl*rho*c*w*u_rel_mag**2)
        drag = (0.5*cd*rho*c*w*u_rel_mag**2)

        # Write the lift and drag magnitude at each actuator node to a CSV file
        write_lift_and_drag('lift', simTime, theta, lift)
        write_lift_and_drag('drag', simTime, theta, drag)

        # Rotate the blade unit vectors to be pointing in the rotated positions
        blade_unit_vec = np.dot(Rz, np.dot(Rx, blade_unit_vec_base))

        # Tile the blade coordinates for every mesh point, [numGridPts*ndim x problem.num_blade_segments]
        blade_pos_full = np.tile(blade_pos, (np.shape(coords)[0], 1))

        # Subtract and square to get the dx^2 values in the x, y, and z directions
        dx_full = (coordsLinear - blade_pos_full)**2

        # Add together to get |x^2 + y^2 + z^2|^2
        dist2 = dx_full[0::ndim] + dx_full[1::ndim] + dx_full[2::ndim]

        # Calculate the force magnitude at every mesh point due to every node [numGridPts x NumActuators]
        nodal_lift = lift*np.exp(-dist2/eps**2)/(eps**3 * np.pi**1.5)
        nodal_drag = drag*np.exp(-dist2/eps**2)/(eps**3 * np.pi**1.5)
            
        for k in range(problem.num_blade_segments):
            # The drag unit simply points opposite the relative velocity unit vector
            drag_unit_vec = -np.copy(u_unit_vec[:, k])
            
            # The lift is normal to the plane generated by the blade and relative velocity
            lift_unit_vec = np.cross(drag_unit_vec, blade_unit_vec[:, 1])

            # lift_unit_mag = np.linalg.norm(lift_unit_vec)
            # if lift_unit_mag < 1e-6:
            #     lift_unit_mag = 1e-6
            # lift_unit_vec = lift_unit_vec/lift_unit_mag

            # All force magnitudes get multiplied by the correctly-oriented unit vector
            vector_nodal_lift = np.outer(nodal_lift[:, k], lift_unit_vec)
            vector_nodal_drag = np.outer(nodal_drag[:, k], drag_unit_vec)

            if dfd == None:
                lift_force += vector_nodal_lift
                drag_force += vector_nodal_drag

            elif dfd == 'c_lift':
                for j in range(ndim):
                    dfd_c_lift[j::ndim, k] += vector_nodal_lift[:, j]

            elif dfd == 'c_drag':
                for j in range(ndim):
                    dfd_c_drag[j::ndim, k] += vector_nodal_drag[:, j]

            elif dfd == 'chord':
                for j in range(ndim):
                    dfd_chord[j::ndim, k] += vector_nodal_lift[:, j] + vector_nodal_drag[:, j]

            # Compute the total force vector [x, y, z] at a single actuator node
            actuator_lift = lift[k]*lift_unit_vec
            actuator_drag = drag[k]*drag_unit_vec

            # Note: since this will be used to define the force (torque) from fluid -> blade
            # we reverse the direction that otherwise gives the turbine force from blade -> fluid
            actuator_force = -(actuator_lift + actuator_drag)

            # Find the component in the direction tangential to the blade
            tangential_actuator_force = np.dot(actuator_force, blade_unit_vec[:, 2])

            # Multiply by the distance away from the hub to get a torque
            actuator_torque = tangential_actuator_force*rdim[k]

            # Add to the total torque
            problem.rotor_torque += actuator_torque

            if save_debugging_files:
                idx = problem.num_blade_segments*blade_ct + k
                blade_pos_paraview[idx, :] = blade_pos[:, k]
                lift_paraview[idx, :] = actuator_lift
                drag_paraview[idx, :] = actuator_drag
                torque_paraview[idx, :] = tangential_actuator_force*blade_unit_vec[:, 2]

    if save_debugging_files:
        fp = open('./output/%s/timeSeries/actuatorForces%05d.vtk' % (problem.params.name, problem.num_times_called), 'w')
        fp.write('# vtk DataFile Version 2.0\n')
        fp.write('Lift, Drag, and Torque Forces\n')
        fp.write('ASCII\n')
        fp.write('DATASET UNSTRUCTURED_GRID\n')

        def write_paraview_vector(name, data):
            num_points = np.shape(data)[0]

            if name is None:
                fp.write('POINTS %d float\n' % (num_points))
            else:        
                fp.write('VECTORS %s float\n' % (name))
            for k in range(num_points):
                x_val = data[k, 0]
                y_val = data[k, 1]
                z_val = data[k, 2]
                if k < num_points-1:
                    fp.write('%.3e %.3e %.3e  ' % (x_val, y_val, z_val))
                else:
                    fp.write('%.3e %.3e %.3e\n\n' % (x_val, y_val, z_val))
                    
            if name is None:
                fp.write('POINT_DATA %d\n' % (num_points))

        write_paraview_vector(None, blade_pos_paraview)
        write_paraview_vector('Lift', lift_paraview)
        write_paraview_vector('Drag', drag_paraview)
        write_paraview_vector('Torque', torque_paraview)

        problem.num_times_called += 1
        fp.close()


    if dfd == None:
        # The total turbine force is the sum of lift and drag effects
        turbine_force = drag_force + lift_force

        # Riffle-shuffle the x-, y-, and z-column force components
        for k in range(ndim):
            tf_vec[k::ndim] = turbine_force[:, k]

        # Remove near-zero values
        tf_vec[np.abs(tf_vec) < 1e-12] = 0.0

        # Initialize turbine force (tf) dolfin function
        tf = Function(problem.fs.V)

        # Set the dolfin vector values
        tf.vector()[:] = tf_vec

        return tf

        # else:
        #     # FIXME: We may need to use .assign() here instead of this
        # FIXED, we assign it in solver manager now
        #     # tf.assign()
        #     tf.vector()[:] = tf_vec

    elif dfd == 'c_lift':
        save_c_lift = False

        if save_c_lift:
            save_derivative_file('dfdcl', dfd_c_lift)

        return dfd_c_lift

    elif dfd == 'c_drag':
        save_c_drag = False

        if save_c_drag:
            save_derivative_file('dfdcd', dfd_c_drag)

        return dfd_c_drag

    elif dfd == 'chord':
        save_chord = True

        if save_chord:
            save_derivative_file('dfdchord', dfd_chord)

        return dfd_chord

#================================================================

def CalculateActuatorLineTurbineForces(problem, simTime, dfd=None, verbose=False):
    if problem.params.get("optimization",{}) and verbose:
        print("Current Optimization Time: "+repr(simTime))

    # if dfd is None, alm_output is a dolfin function (tf) [1 x numPts*ndim]
    # otherwise, it returns a numpy array of derivatives [numPts*ndim x numControls]
    alm_output = UpdateActuatorLineForce(problem, simTime, dfd)

    return alm_output

#================================================================

