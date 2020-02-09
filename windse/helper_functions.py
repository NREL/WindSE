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
            temp[dim*(sparse_ids)+0] = np.sum(actuators_x*n[i],axis=1)
            temp[dim*(sparse_ids)+1] = np.sum(actuators_y*n[i],axis=1)
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
            tfs[i][0::dim] = d_actuators_x*n[i]
            tfs[i][1::dim] = d_actuators_y*n[i]

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
            tfs[i][0::dim] = d_actuators_x*n[i]
            tfs[i][1::dim] = d_actuators_y*n[i]

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
            tfs[i][0::dim] = d_actuators_x*n[i]+actuators_x*d_n[i]
            tfs[i][1::dim] = d_actuators_y*n[i]+actuators_y*d_n[i]

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
