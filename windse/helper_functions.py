import __main__
import os, sys

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from windse import windse_parameters
    if windse_parameters.dolfin_adjoint:
        from dolfin import dx, File, dot
        from dolfin_adjoint import Constant, Function, Expression, assemble
    else:
        from dolfin import Constant, Function, Expression, dot, dx, assemble, File

    import numpy as np
    import scipy.interpolate as interp
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

def UpdateActuatorLineForce(problem, u_local, simTime_id, dt, turb_i, dfd=None, verbose=False):

    simTime = problem.simTime_list[simTime_id]


    if verbose:
        print("Current Optimization Time: "+repr(simTime)+", Turbine #"+repr(turb_i))
        sys.stdout.flush()

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

    def save_derivative_file(folder,filename, deriv_array):

        dolfin_function = Function(problem.fs.V)

        fp = File(folder+'%s.pvd' % (filename))

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


    def build_lift_and_drag(problem, u_rel, blade_unit_vec, rdim, twist, c):

        def get_angle_between_vectors(a, b, n):
            a_x_b = np.dot(np.cross(n, a), b)

            norm_a = np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
            norm_b = np.sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])

            c1 = a_x_b/(norm_a*norm_b)
            c1 = np.clip(c1, -1.0, 1.0)
            aoa_1 = np.arcsin(c1)

            c2 = np.dot(a, b)/(norm_a*norm_b)
            c2 = np.clip(c2, -1.0, 1.0)
            aoa_2 = np.arccos(c2)
            
            if aoa_2 > np.pi/2.0:
                if aoa_1 < 0:
                    aoa_1 = -np.pi - aoa_1
                else:
                    aoa_1 = np.pi - aoa_1
            
            aoa_1_deg = aoa_1/np.pi*180.0
            
            return aoa_1


        # If this is the first time calling the function...
        if problem.first_call_to_alm:
            # build the lift-drag table interpolators
            rdim_all = np.linspace(0, rdim[-1], np.shape(problem.lift_table)[1])
            problem.interp_lift = interp.RectBivariateSpline(problem.interp_angles, rdim_all, problem.lift_table)
            problem.interp_drag = interp.RectBivariateSpline(problem.interp_angles, rdim_all, problem.drag_table)


        # Initialize the real cl and cd profiles
        real_cl = np.zeros(problem.num_blade_segments)
        real_cd = np.zeros(problem.num_blade_segments)

        fp = open(problem.aoa_file, 'a')

        tip_loss = np.zeros(problem.num_blade_segments)

        for k in range(problem.num_blade_segments):
            # Get the relative wind velocity at this node
            wind_vec = u_rel[:, k]

            # Remove the component in the radial direction (along the blade span)
            wind_vec -= np.dot(wind_vec, blade_unit_vec[:, 1])*blade_unit_vec[:, 1]

            # aoa = get_angle_between_vectors(arg1, arg2, arg3)
            # arg1 = in-plane vector pointing opposite rotation (blade sweep direction)
            # arg2 = relative wind vector at node k, including blade rotation effects (wind direction)
            # arg3 = unit vector normal to plane of rotation, in this case, radially along span
            aoa = get_angle_between_vectors(-blade_unit_vec[:, 2], wind_vec, -blade_unit_vec[:, 1])

            # Compute tip-loss factor
            if rdim[k] < 1e-12:
                tip_loss[k] = 1.0
            else:
                loss_exponent = 3.0/2.0*(rdim[-1]-rdim[k])/(rdim[k]*np.sin(aoa))
                acos_arg = np.exp(-loss_exponent)
                acos_arg = np.clip(acos_arg, -1.0, 1.0)
                tip_loss[k] = 2.0/np.pi*np.arccos(acos_arg)

            # Remove the portion of the angle due to twist
            aoa -= twist[k]


            # if np.abs(problem.farm.baseline_chord[k] - c[k]) > 1e-6:
            #     # print('Chord differs from baseline at node %d:' % (k))
            #     c_factor = c[k]/problem.farm.baseline_chord[k]

            #     # Penalize by 10 degrees if the chord is doubled
            #     penalty = c_factor*(np.pi*5.0/180.0)

            #     aoa -= penalty

            #     # real_cd[k] *= c_factor


            # Store the cl and cd by interpolating this (aoa, span) pair from the tables
            real_cl[k] = problem.interp_lift(aoa, rdim[k])
            real_cd[k] = problem.interp_drag(aoa, rdim[k])



            # Write the aoa to a file for future reference
            fp.write('%.5f, ' % (aoa/np.pi*180.0))

        fp.close()

        return real_cl, real_cd, tip_loss


    #================================================================
    # Get Mesh Properties
    #================================================================

    ndim = problem.dom.dim

    # Get the coordinates of the vector function space
    coords = problem.fs.V.tabulate_dof_coordinates()
    coords = np.copy(coords[0::problem.dom.dim, :])

    # Resape a linear copy of the coordinates for every mesh point
    coordsLinear = np.copy(coords.reshape(-1, 1))

    # Initialize a single turbine force function (contains only a single turbine for power calculations)
    tf_individual = Function(problem.fs.V)

    # Initialize a cumulative turbine force function (contains all turbines)
    tf = Function(problem.fs.V)
    tf.vector()[:] = 0.0

    gaussian_fn = Function(problem.fs.Q)
    gaussian_fn.vector()[:] = 0.0

    # Initialize a cylindrical field function
    cyld = Function(problem.fs.V)

    #================================================================
    # Set Turbine and Fluid Properties
    #================================================================

    # Set the density
    rho = 1.0

    # Set the number of blades in the turbine
    num_blades = 3

    # Width of Gaussian
    # Note: this sets the gaussian width to roughly twice the minimum cell length scale
    eps = problem.gaussian_width

    # Initialize the sum of the torque on the rotor shaft to zero
    # problem.rotor_torque = 0.0
    if not hasattr(problem,"rotor_torque"):
        problem.rotor_torque = np.zeros(problem.farm.numturbs)
        problem.rotor_torque_dolfin = np.zeros(problem.farm.numturbs)


    # initialize numpy torque
    rotor_torque_numpy_temp = 0.0

    # for turb_i in range(problem.farm.numturbs):
    # for turb_i in turb_index:
    # print(turb_i)

    # # Set the hub height
    # hub_height = problem.farm.HH[turb_i]

    # # Set the yaw of the turbine
    # yaw = problem.farm.yaw[turb_i]

    # Blade length (turbine radius)
    L = problem.farm.radius[turb_i]

    #================================================================
    # Set Derived Constants
    #================================================================

    # Calculate the radial position of each actuator node
    rdim = np.linspace(0.0, L, problem.num_blade_segments)

    # Calculate width of an individual blade segment
    # w = rdim[1] - rdim[0]
    w = (rdim[1] - rdim[0])*np.ones(problem.num_blade_segments)
    w[0] = w[0]/2.0
    w[-1] = w[-1]/2.0

    # Calculate an array describing the x, y, z position of each actuator node
    # Note: The basic blade is oriented along the +y-axis
    blade_pos_base = np.vstack((np.zeros(problem.num_blade_segments),
                        rdim,
                        np.zeros(problem.num_blade_segments)))

    # Calculate the blade velocity
    angular_velocity = 2.0*np.pi*problem.rpm/60.0
    tip_speed = angular_velocity*L
    # tip_speed = 9.0

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
    # blade_unit_vec_base[:, 2] = points tangential to blade span axis (generates a torque about rotor shaft)
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
    twist = np.array(problem.mtwist[turb_i], dtype = float)

    cl = np.array(problem.mcl[turb_i], dtype = float)
    cd = np.array(problem.mcd[turb_i], dtype = float)
    tip_loss = np.ones(problem.num_blade_segments)


    # if problem.first_call_to_alm:
    #     print('twist', twist)
    #     print('lift-table', problem.lift_table)
    #     print('drag-table', problem.drag_table)

    # cl = np.ones(problem.num_blade_segments)
    # cd = np.ones(problem.num_blade_segments)

    # Read the chord length from the values specified in the problem manager
    # c = L/20.0
    c = np.array(problem.mchord[turb_i], dtype = float)


    # print()
    # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # print("current controls from helper")
    # print(np.array(c,dtype=float))
    # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # print()


    # Initialze arrays depending on what this function will be returning
    if dfd is None:
        tf_vec = np.zeros(np.size(coords))
        tf_vec_for_power = np.zeros(np.size(coords))
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

    save_safety_switch_files = False

    if save_safety_switch_files:
        blade_pos_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        lift_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        drag_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        torque_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        torque_force_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        u_fluid_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        u_blade_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        u_rel_paraview = np.zeros((problem.num_blade_segments*num_blades, 3))
        chord_paraview = np.hstack((c, c, c))

        # print('k-1/2, measuring fluid at  : ', 0.5*(problem.simTime_list[-1] + simTime))
        # print('k    , simTime:              ', simTime)
        # print('k+1/2, calculating force at: ', simTime+0.5*problem.dt)
        # print('k+1  , next sim Time:        ', simTime+problem.dt)


    # Calculate the blade position based on current simTime and turbine RPM
    period = 60.0/problem.rpm
    # theta_offset = simTime/period*2.0*np.pi
    theta_offset = (simTime+0.5*dt)/period*2.0*np.pi
    # theta_offset = 0.0

    baseline_chord = problem.farm.baseline_chord


    baseleline_chord_vols = np.pi*(baseline_chord/2.0)**2*w
    c_vols = np.pi*(c/2.0)**2*w
    # baseleline_chord_vols = 4.0/3.0*np.pi*(baseline_chord/2.0)**3
    # c_vols = 4.0/3.0*np.pi*(c/2.0)**3

    c0 = np.dot(baseleline_chord_vols, rdim)
    c1 = np.dot(c_vols, rdim)
    baseline_omega = 2.0*np.pi*10.6/60.0
    w_new = np.sqrt((baseline_omega**2*c0)/c1)
    rpm_new = w_new*60.0/(2.0*np.pi)

    # print(c)


    # print('Turbine currently spinning at %.3f RPM' % (problem.rpm))
    # print('Turbine should be spinning at %.3f RPM' % (rpm_new))
    # print('To match centripetal forces from baseline case.')

    print("Current Yaw: "+repr(float(problem.farm.myaw[turb_i])))

    # Treat each blade separately
    for blade_ct, theta_0 in enumerate(theta_vec):
        theta = theta_0 + theta_offset

        # Generate a rotation matrix for this turbine blade
        Rx = rot_x(theta)
        Rz = rot_z(float(problem.farm.myaw[turb_i]))

        # Rotate the blade velocity in the global x, y, z, coordinate system
        # Note: blade_vel_base is negative since we seek the velocity of the fluid relative to a stationary blade
        # and blade_vel_base is defined based on the movement of the blade
        blade_vel = np.dot(Rz, np.dot(Rx, -blade_vel_base))

        # Rotate the blade unit vectors to be pointing in the rotated positions
        blade_unit_vec = np.dot(Rz, np.dot(Rx, blade_unit_vec_base))

        # Rotate the entire [x; y; z] matrix using this matrix, then shift to the hub location
        blade_pos = np.dot(Rz, np.dot(Rx, blade_pos_base))
        blade_pos[0, :] += problem.farm.x[turb_i]
        blade_pos[1, :] += problem.farm.y[turb_i]
        blade_pos[2, :] += problem.farm.z[turb_i]

        # time_offset = 1
        # if len(problem.simTime_list) < time_offset:
        #     theta_behind = theta_0 + problem.simTime_list[0]/period*2.0*np.pi
        # else:
        #     theta_behind = theta_0 + problem.simTime_list[-time_offset]/period*2.0*np.pi

        problem.first_call_to_function = False
        time_offset = 1
        if simTime_id < time_offset:
            theta_behind = theta_0 + 0.5*(problem.simTime_list[simTime_id]+simTime)/period*2.0*np.pi
        else:
            theta_behind = theta_0 + 0.5*(problem.simTime_list[simTime_id-time_offset]+simTime)/period*2.0*np.pi
            # if blade_ct == 0:
            #     print('SimTime = %f, using %f' % (simTime, problem.simTime_list[-time_offset]))

        Rx_alt = rot_x(theta_behind)
        blade_pos_alt = np.dot(Rz, np.dot(Rx_alt, blade_pos_base))
        blade_pos_alt[0, :] += problem.farm.x[turb_i]
        blade_pos_alt[1, :] += problem.farm.y[turb_i]
        blade_pos_alt[2, :] += problem.farm.z[turb_i]

        # print('turbine %d x-pos = %f, y-pos = %f, hh = %f, yaw = %f' % (turb_i,
        #     problem.farm.x[turb_i], problem.farm.y[turb_i], problem.farm.HH[turb_i], problem.farm.yaw[turb_i]))

        # Initialize space to hold the fluid velocity at each actuator node
        u_fluid = np.zeros((3, problem.num_blade_segments))

        # Set if using local velocity around inidividual nodes
        using_local_velocity = True

        if using_local_velocity:
            # Generate the fluid velocity from the actual node locations in the flow
            for k in range(problem.num_blade_segments):

                u_fluid[:, k] = u_local(blade_pos_alt[0, k],
                         blade_pos_alt[1, k],
                         blade_pos_alt[2, k])

                u_fluid[:, k] -= np.dot(u_fluid[:, k], blade_unit_vec[:, 1])*blade_unit_vec[:, 1]

        else:
            # Generate the fluid velocity using the inlet velocity (x-pos = x_range[0])        
            for k in range(problem.num_blade_segments):
                # u_fluid[:, k] = [problem.bd.HH_vel,0,0]

                u_fluid[:, k] = problem.bd.bc_velocity(problem.dom.x_range[0],
                                             blade_pos_alt[1, k],
                                             blade_pos_alt[2, k])

                u_fluid[1, k] = 0.0
                u_fluid[2, k] = 0.0
                # Force inlet velocity (1, 0, 0) for safe_mode
                # u_fluid[:, k] = np.array([1, 0, 0])

        # print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        # print(u_fluid)
        # print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

        problem.blade_pos_previous[blade_ct] = blade_pos

                        
        # Form the total relative velocity vector (including velocity from rotating blade)
        u_rel = u_fluid + blade_vel

        # Force fluid-only relative velocity for safe_mode
        # u_rel = u_fluid

        # Create unit vectors in the direction of u_rel


        # u_rel_mag = np.sqrt(u_rel[0]**2.0+u_rel[1]**2.0+u_rel[2]**2.0)
        # u_rel_mag = u_rel[0]**2.0+u_rel[1]**2.0+u_rel[2]**2.0

        u_rel_mag = np.linalg.norm(u_rel, axis=0)
        u_rel_mag[u_rel_mag < 1e-6] = 1e-6
        u_unit_vec = u_rel/u_rel_mag
        
        # Calculate the lift and drag forces using the relative velocity magnitude
        # if dfd == "u":
        #     d_u_rel_mag = 2*u_rel[0]+2*u_rel[1]+2*u_rel[2]
        #     lift = (0.5*cl*rho*c*w*2*u_rel_mag*d_u_rel_mag)
        # else:
        using_lift_and_drag_tables = True

        if using_lift_and_drag_tables:
            cl, cd, tip_loss = build_lift_and_drag(problem, u_rel, blade_unit_vec, rdim, twist, c)

        
        # Calculate the lift and drag forces using the relative velocity magnitude
        lift = tip_loss*(0.5*cl*rho*c*w*u_rel_mag**2)
        drag = tip_loss*(0.5*cd*rho*c*w*u_rel_mag**2)

        # Write the lift and drag magnitude at each actuator node to a CSV file
        write_lift_and_drag('lift', simTime, theta, lift)
        write_lift_and_drag('drag', simTime, theta, drag)


        # Tile the blade coordinates for every mesh point, [numGridPts*ndim x problem.num_blade_segments]
        blade_pos_full = np.tile(blade_pos, (np.shape(coords)[0], 1))

        # Subtract and square to get the dx^2 values in the x, y, and z directions
        dx_full = (coordsLinear - blade_pos_full)**2

        # Add together to get |x^2 + y^2 + z^2|^2
        dist2 = dx_full[0::ndim] + dx_full[1::ndim] + dx_full[2::ndim]


        gaussian_field = np.exp(-dist2/eps**2)/(eps**3 * np.pi**1.5)
        gauss_add = np.sum(gaussian_field, axis=1)
        gauss_add[np.abs(gauss_add) < 1e-12] = 0.0
        gaussian_fn.vector()[:] += gauss_add

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

            # print('Scalar Lift/Summed Lift = %.6e' % (lift[k]/np.sum(nodal_lift[:, k])))
            # print('Scalar Drag/Summed Drag = %.6e' % (drag[k]/np.sum(nodal_drag[:, k])))
            # print('sum_gauss = %.6e' % (np.sum(gaussian_field[:, k])))

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
            # actuator_force = -(actuator_lift - actuator_drag)

            # Find the component in the direction tangential to the blade
            tangential_actuator_force = np.dot(actuator_force, blade_unit_vec[:, 2])

            # Multiply by the distance away from the hub to get a torque
            actuator_torque = tangential_actuator_force*rdim[k]

            # Add to the total torque
            rotor_torque_numpy_temp += actuator_torque  ### Should this be an output?

            if save_safety_switch_files:
                idx = problem.num_blade_segments*blade_ct + k
                blade_pos_paraview[idx, :] = blade_pos[:, k]
                lift_paraview[idx, :] = actuator_lift
                drag_paraview[idx, :] = actuator_drag
                torque_paraview[idx, :] = tangential_actuator_force*blade_unit_vec[:, 2]
                torque_force_paraview[idx, :] = -actuator_force

                u_rel_paraview[idx, :] = u_rel[:, k]
                u_blade_paraview[idx, :] = blade_vel[:, k]
                u_fluid_paraview[idx, :] = u_fluid[:, k]

    # Output the numpy version of rotor_torque
    problem.rotor_torque[turb_i] = rotor_torque_numpy_temp


    if save_safety_switch_files:
        folder_string = problem.params.folder+"timeSeries/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)
        fp = open(folder_string +('actuatorForces%05d.vtk' % (problem.num_times_called)), 'w')
        fp.write('# vtk DataFile Version 2.0\n')
        fp.write('Lift, Drag, and Torque Forces\n')
        fp.write('ASCII\n')
        fp.write('DATASET UNSTRUCTURED_GRID\n')

        def write_paraview_vector(name, data, scalar_data=False):
            num_points = np.shape(data)[0]

            if scalar_data:
                fp.write('SCALARS %s float\n' % (name))
                fp.write('LOOKUP_TABLE default\n')

                for k in range(num_points):
                    if k < num_points-1:
                        fp.write('%.3e ' % (data[k]))
                    else:
                        fp.write('%.3e \n\n' % (data[k]))

            else:
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
        write_paraview_vector('Chord', chord_paraview, scalar_data=True)

        write_paraview_vector('Lift', lift_paraview)
        write_paraview_vector('Drag', drag_paraview)
        write_paraview_vector('Torque', torque_paraview)
        write_paraview_vector('Torque_force', torque_force_paraview)

        write_paraview_vector('u_rel', u_rel_paraview)
        write_paraview_vector('u_fluid', u_fluid_paraview)
        write_paraview_vector('u_blade', u_blade_paraview)

        fp.close()

    if dfd == None:
        # The total turbine force is the sum of lift and drag effects
        turbine_force = drag_force + lift_force
        turbine_force_for_power = -drag_force + lift_force

        # Riffle-shuffle the x-, y-, and z-column force components
        for k in range(ndim):
            tf_vec[k::ndim] = turbine_force[:, k]
            tf_vec_for_power[k::ndim] = turbine_force_for_power[:, k]

        # Remove near-zero values
        tf_vec[np.abs(tf_vec) < 1e-12] = 0.0

        # Assign the individual turbine force
        tf_individual.vector()[:] = tf_vec#_for_power

        # Create a cylindrical expression aligned with the position of this turbine
        cyld_expr = Expression(('sin(yaw)*(x[2]-zs)', '-cos(yaw)*(x[2]-zs)', '(x[1]-ys)*cos(yaw)-(x[0]-xs)*sin(yaw)'),
            degree=1,
            yaw=problem.farm.myaw[turb_i],
            xs=problem.farm.mx[turb_i],
            ys=problem.farm.my[turb_i],
            zs=problem.farm.z[turb_i])

        cyld.interpolate(cyld_expr)

        # Integrate dot(Force, In-Plane-Vector) to get total torque
        problem.rotor_torque_dolfin[turb_i] = assemble(dot(-tf_individual, cyld_expr)*dx)

        # Add to the cumulative turbine force
        tf.vector()[:] += tf_vec

        # if turb_i == 0:
        #     # Initialize turbine force (tf) dolfin function
        #     tf = Function(problem.fs.V)

        #     # Set the dolfin vector values
        #     tf.vector()[:] = tf_vec
            
        # else:
        #     tf.vector()[:] += tf_vec

        # Test storing them here for objective function calc
        # problem.tf_individual = tf_individual
        problem.cyld_expr_list[turb_i] = cyld_expr
        problem.cyld = cyld
        # print('just saved turb %.0d at x = %f, y = %.0f, yaw = %.4f' % (turb_i, problem.farm.x[turb_i], problem.farm.y[turb_i], problem.farm.yaw[turb_i]))

    fp = open(problem.aoa_file, 'a')
    fp.write('\n')
    fp.close()

    problem.num_times_called += 1
    problem.first_call_to_alm = False

    # fp = File('%s/timeSeries/gauss.pvd' % (problem.params.folder))
    # gaussian_fn.rename('gauss', 'gauss')
    # fp << (gaussian_fn, len(problem.simTime_list))

    if dfd == None:

        return tf

    elif dfd == 'c_lift':
        save_c_lift = False

        if save_c_lift:
            save_derivative_file(problem.params.folder+"timeSeries/",'dfdcl', dfd_c_lift)

        return dfd_c_lift

    elif dfd == 'c_drag':
        save_c_drag = False

        if save_c_drag:
            save_derivative_file(problem.params.folder+"timeSeries/",'dfdcd', dfd_c_drag)

        return dfd_c_drag

    elif dfd == 'chord':
        save_chord = False

        if save_chord:
            save_derivative_file(problem.params.folder+"timeSeries/",'dfdchord', dfd_chord)

        return dfd_chord

