import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from windse import windse_parameters
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin import dx, File
        from dolfin_adjoint import Constant, Function, assemble
    else:
        from dolfin import Constant, Function, dx, assemble

    import numpy as np

def BaseHeight(x,y,ground,dx=0,dy=0):
    return Constant(ground(float(x),float(y),dx=dx,dy=dy))

def SineProfile(xrot,a,R,dfd=None,d_xrot=None):
    ''' Returns a numpy vector corresponding to a scaler function space'''

    ### Collect the relevant turbine data ###
    A = np.pi*R**2.0

    ### Precompute functions and derivatives ###
    r    = np.sqrt(SquaredRadius_yz(xrot))
    ct   = 4.*0.5*A*a/(1.-a)
    sine = (r/R*np.sin(np.pi*r/R)+0.5)/(.81831)
    if d_xrot is not None:
        d_r    = 0.5*SquaredRadius_yz(xrot,d_xrot=d_xrot)/r
        d_sine = (r/R*np.cos(np.pi*r/R)*(np.pi*d_r/R)+d_r/R*np.sin(np.pi*r/R))/(.81831)
    if dfd == "a":
        d_ct   = 4*0.5*A/(a-1.)**2.0

    ### Calculate output
    if dfd is None:
        out = ct*sine
    elif dfd == "a":
        out = d_ct*sine
    elif d_xrot is not None:
        out = ct*d_sine
    else:
        out = np.zeros(x.shape)

    return out

def ConstantProfile(xrot,a,R,dfd=None,d_xrot=None):
    ''' Returns a numpy vector corresponding to a scaler function space'''

    ### Collect the relevant turbine data ###
    A = np.pi*R**2.0

    ### Calculate output
    if dfd is None:
        out = 4*0.5*A*a/(1.-a)
    elif dfd == "a":
        out = 4*0.5*A/(a-1.)**2.0
    else:
        out = np.zeros(xrot[0].shape)

    return out

def SquaredRadius_yz(xrot,d_xrot=None):
    ### Compute the output and derivatives using the chain rule ###
    if d_xrot is None:
        out = np.power(xrot[1],2.0)+np.power(xrot[2],2.0)
    else:
        out = 2.0*xrot[1]*d_xrot[1]+2.0*xrot[2]*d_xrot[2]
    return out 

def ThicknessKernal(xrot,W,d_xrot=None):
    ### Define integration constant ###
    T_norm = 1.855438667500383

    ### Precompute functions and derivatives ###
    f = -np.power((xrot[0]/W),6.0)
    if d_xrot is not None:
        d_f = -6.0*np.power((xrot[0]/W),5.0) * (d_xrot[0]/W)

    ### Compute the output and derivatives using the chain rule ###
    if d_xrot is None:
        out = np.exp(f)/(T_norm*W)
    else:
        out = d_f * np.exp(f)/(T_norm*W)

    return out

def CylinderKernal(xrot,R,d_xrot=None):
    ### Define integration constant ###
    C_norm = 2.914516237206873

    ### Precompute functions and derivatives ###
    r = SquaredRadius_yz(xrot)
    f = -np.power(r/R**2.0,6.0)
    if d_xrot is not None:
        d_r = SquaredRadius_yz(xrot,d_xrot=d_xrot)
        d_f = -6.0*np.power(r/R**2.0,5.0) * (d_r/R**2.0)

    ### Compute the output and derivatives using the chain rule ###
    if d_xrot is None:
        out = np.exp(f)/(C_norm*R**2.0)
    else:
        out = d_f * np.exp(f)/(C_norm*R**2.0)

    return out

def DiskKernal(xrot, W, R, d_xrot=None):
    ''' Returns a numpy vector corresponding to a scaler function space'''

    ### Precompute functions and derivatives ###
    T = ThicknessKernal(xrot,W)
    C = CylinderKernal(xrot,R)
    if d_xrot is not None:
        d_T = ThicknessKernal(xrot,W,d_xrot=d_xrot)
        d_C = CylinderKernal(xrot,R,d_xrot=d_xrot)

    ### Compute the output and derivatives using the product rule ###
    if d_xrot is None:
        out = T*C
        # out = C
    else:
        out = T*d_C + d_T*C
        # out = d_C

    return out

def Transform(x,x0,HH,yaw,ground,dfd=None):

    ### Get the dimension ###
    dim = x.shape[0]

    ### Rotate and Shift the kernel ###
    if dfd is None:
        xrot =   np.cos(yaw)*(x[0]-x0[0]) + np.sin(yaw)*(x[1]-x0[1])
        yrot = - np.sin(yaw)*(x[0]-x0[0]) + np.cos(yaw)*(x[1]-x0[1])
        if dim == 3:
            zrot =   x[2]-(ground(x0[0],x0[1])+HH)
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
        xrot = - np.sin(yaw)*(x[0]-x0[0]) + np.cos(yaw)*(x[1]-x0[1])
        yrot = - np.cos(yaw)*(x[0]-x0[0]) - np.sin(yaw)*(x[1]-x0[1])
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

def ActuatorDisk(x,wind_farm,i,fs,dfd=None):
    ''' Returns a numpy vector corresponding to a vector function space '''

    ### Pick the correct force function ###
    if wind_farm.force == "constant":
        ForceFunction = ConstantProfile
    elif wind_farm.force == "sine":
        ForceFunction = SineProfile
    else:
        ValueError("Unknown force type: "+wind_farm.force)

    ### Collect the relevant turbine data ###
    x0 = [float(wind_farm.mx[i]),float(wind_farm.my[i]),float(wind_farm.mz[i])]
    yaw = float(wind_farm.myaw[i])+wind_farm.inflow_angle
    a = float(wind_farm.ma[i])
    HH = wind_farm.HH[i]
    W = wind_farm.W[i]
    R = wind_farm.RD[i]/2.0
    dim = x.shape[0]

    # ### Create Empty Functions ###
    # actuator = Function(fs.tf_V)
    # actuator_x = Function(fs.tf_V0)
    # actuator_y = Function(fs.tf_V1)
    actuator = np.zeros(x.shape[1]*3)

    ### Precompute functions and derivatives ###
    xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground)
    force = ForceFunction(xrot,a,R)
    space_kernal = DiskKernal(xrot,W,R)
    disk = force*space_kernal
    # disk = space_kernal
    # disk = x0[0]**2
    # disk = np.exp(-np.power(xrot[0],2.0))
    # W = 0.02
    # disk = np.exp(-np.power(xrot[0]/W,6.0))
    # disk = np.exp(-xrot[1])
    if dfd is not None:
        d_xrot = Transform(x,x0,HH,yaw,wind_farm.dom.Ground,dfd=dfd)
        d_force = ForceFunction(xrot,a,R,dfd=dfd,d_xrot=d_xrot)
        d_space_kernal = DiskKernal(xrot,W,R,d_xrot=d_xrot)
        d_disk = force*d_space_kernal + d_force*space_kernal
        # d_disk = d_space_kernal
        # d_disk = 2.0*xrot[0]*d_xrot[0]
        # d_disk = -2.0*np.exp(-np.power(xrot[0],2.0))*(xrot[0]*d_xrot[0])
        # d_disk = -np.exp(-np.power(xrot[0]/W,6.0))*6.0*(np.power(xrot[0]/W,5.0))*(d_xrot[0]/W)#+2.0*xrot[0]*d_xrot[0]*np.exp(-np.power(xrot[0]/W,6.0))
        # if dfd == "x":
        #     d_disk = 2.0*x0[0]
        # else:
        #     d_disk = 0.0


    ### Create the normal function and it's derivatives ###
    nx = np.cos(yaw)
    ny = np.sin(yaw)
    if dfd == "yaw":
        d_nx = -1*np.sin(yaw)
        d_ny = np.cos(yaw)
    else:
        d_nx = 0.0
        d_ny = 0.0

    ### Compute the output and derivatives using the product rule ###
    if dfd is None:
        # actuator_x.vector()[:] = disk*nx
        # actuator_y.vector()[:] = disk*ny
        actuator[0::dim] = disk*nx
        actuator[1::dim] = disk*ny
    else:
        # actuator_x.vector()[:] = d_disk*nx + disk*d_nx
        # actuator_y.vector()[:] = d_disk*ny + disk*d_ny
        actuator[0::dim] = d_disk*nx + disk*d_nx
        actuator[1::dim] = d_disk*ny + disk*d_ny

    ### Convert scaler to vector ####
    # if dim == 3:
    #     actuator_z = Function(fs.tf_V2)
    #     fs.tf_Assigner.assign(actuator,[actuator_x,actuator_y,actuator_z])
    # else:
    #     fs.tf_Assigner.assign(actuator,[actuator_x,actuator_y])

    ### Filter out values smaller than machine precision ###
    # filtered = actuator.vector().get_local()
    tol = 1e-17
    # filtered[abs(filtered) < tol] = 0.0
    actuator[abs(actuator) < tol] = 0.0



    # if dfd is None:
    #     actuator[0::dim] = xrot[0]
    #     actuator[1::dim] = xrot[0]
    #     actuator[2::dim] = xrot[0]
    # else:
    #     actuator[0::dim] = d_xrot[0]
    #     actuator[1::dim] = d_xrot[0]
    #     actuator[2::dim] = d_xrot[0]





    ### Return just the numpy array ###
    # return filtered
    return actuator

def SingleDiskTurbineForces(x,wind_farm,i,fs,dfd=None,save_actuator=False):
    ''' Returns a numpy vector corresponding to a vector function space '''

    ### Collect the relevant turbine data ###
    yaw = float(wind_farm.myaw[i])+wind_farm.inflow_angle

    # print()
    # print("inflow angle: " + repr(wind_farm.inflow_angle))
    # print("yaw:          " + repr(float(wind_farm.myaw[i])))
    # print("adjusted yaw: " + repr(yaw))
    # print()


    ### Precompute functions and derivatives ###
    actuator = ActuatorDisk(x,wind_farm,i,fs)
    if dfd is not None:
        d_actuator = ActuatorDisk(x,wind_farm,i,fs,dfd=dfd) 

    ### Calculate the normal squared components and their derivatives ###
    n1 = np.cos(yaw)**2
    n2 = np.sin(yaw)**2
    n3 = 2.0*np.cos(yaw)*np.sin(yaw)
    if dfd == "yaw":
        d_n1 = (-2)*np.cos(yaw)*np.sin(yaw)
        d_n2 = 2*np.sin(yaw)*np.cos(yaw)
        d_n3 = 2.0*(np.cos(2*yaw))
    else:
        d_n1 = 0.0
        d_n2 = 0.0
        d_n3 = 0.0

    ### Compute the output and derivatives using the product rule ###
    if dfd is None:
        tf1 = actuator*n1
        tf2 = actuator*n2
        tf3 = actuator*n3
    else:
        tf1 = d_actuator*n1 + actuator*d_n1    
        tf2 = d_actuator*n2 + actuator*d_n2 
        tf3 = d_actuator*n3 + actuator*d_n3

    # print()
    # print()
    # print(dfd)
    # testf = Function(fs.tf_V0)
    # testf.vector()[:]=1.0
    # test = assemble(testf*dx)
    # print(test)
    # print(sum(tf1[0::3]))
    # print(sum(tf1[0::3])*test)
    # print(yaw)
    # print(max(tf1))
    # print(min(tf1))
    # print()
    # print()


    if save_actuator:
        if dfd is not None:
            return [tf1,tf2,tf3,d_actuator]
        else:
            return [tf1,tf2,tf3,actuator]
    else:
        return [tf1,tf2,tf3]

def CalculateCombinedTurbineForces(x,wind_farm,fs,dfd=None,df_index=None,save_actuator=False):

    print("building TF with: "+repr(dfd))
    ### Initialize the dolfin functions ###
    tf1  = Function(fs.tf_V)
    tf2  = Function(fs.tf_V)
    tf3  = Function(fs.tf_V)
    actuator_disks = []

    if dfd is None:
        ### Calculate all wind_farm and collect them into single functions ###
        for i in range(wind_farm.numturbs):
            out = SingleDiskTurbineForces(x,wind_farm,i,fs,dfd=dfd,save_actuator=save_actuator)
            tf1.vector()[:] += out[0]
            tf2.vector()[:] += out[1]
            tf3.vector()[:] += out[2]
            if save_actuator:
                actuator_disk = Function(fs.tf_V)
                actuator_disk.vector()[:] = out[3]
                actuator_disks.append(actuator_disk)

        # print(dfd)
        # print(assemble(-tf1[0]*dx))
        # File("tf1.pvd") << tf1


    else:
        ### Calculate the derivative vectors ###
        out = SingleDiskTurbineForces(x,wind_farm,df_index,fs,dfd=dfd)
        tf1.vector()[:] += out[0]
        tf2.vector()[:] += out[1]
        tf3.vector()[:] += out[2]

        # print(dfd)
        # print(assemble(-tf1[0]*dx))
        # File("dtf1_d"+dfd+".pvd") << tf1


    if save_actuator:
        return [tf1,tf2,tf3,actuator_disks]
    else:
        return [tf1,tf2,tf3]
