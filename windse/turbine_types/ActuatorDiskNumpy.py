# import windse function
from . import ActuatorDisk
from windse import windse_parameters
if windse_parameters.dolfin_adjoint:
    from windse.blocks import blockify, ActuatorDiskForceBlock

# import dolfin functions
from . import cos, sin, exp, sqrt, dot, Function

import numpy as np
from scipy.special import gamma

class ActuatorDiskNumpy(ActuatorDisk):

    def __init__(self, i,x,y,dom,imported_params=None):
        # initialize the rest of the turbine
        super(ActuatorDiskNumpy, self).__init__(i,x,y,dom,imported_params)

        # blockify custom functions so dolfin adjoint can track them
        if self.params.performing_opt_calc:
            block_kwargs = {
                "construct_sparse_ids":self.construct_sparse_ids,
                "control_types": self.params["optimization"]["control_types"],
                "turb": self
            }
            self.build_actuator_disk = blockify(self.build_actuator_disk,ActuatorDiskForceBlock,block_kwargs=block_kwargs)

    def construct_sparse_ids(self, x, sparse_RDs):
        """ 
        Determines relevant dofs that will be nonzero
        """

        ### Create the bounding box for each turbine ###
        dim = x.shape[0]
        x0 = np.array([[self.mx],[self.my],[self.mz]],dtype=float)
        bounding_limit = [sparse_RDs*self.RD/2.0]*3
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

        return sparse_ids

    ### TODO: this can eventually be replaced by whatever scaffolding we come up with for floating turbines
    def transform(self,x,x0,HH,yaw,ground,dfd=None):

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

    def build_actuator_disk(self,x,inflow_angle,fs,sparse_ids,actuator_disk=None,dfd=None):

        ### Collect the relevant turbine data ###
        x0 = np.array([[self.mx],[self.my],[self.mz]],dtype=float)
        yaw = np.array(self.myaw,dtype=float)+inflow_angle
        a = np.array(self.maxial,dtype=float)
        HH = self.HH
        W = self.thickness
        R = self.RD/2.0

        ### Select sparse x values ###
        dim, N = x.shape
        x=x[:,sparse_ids]

        ### Set up some dim dependent values ###
        S_norm = (2.0+np.pi)/(2.0*np.pi)
        T_norm = 2.0*gamma(7.0/6.0)
        if dim == 3:
            A = np.pi*R**2.0 
            D_norm = np.pi*gamma(4.0/3.0)
        else:
            A = 2*R 
            D_norm = 2.0*gamma(7.0/6.0)
        volNormalization = T_norm*D_norm*W*R**(dim-1)

        ### Define Radial Force Functions ###
        if self.force == "constant":
            def RForce(r): return -1.0
            def dRForce(r,d_r): return 0.0
        elif self.force == "sine":
            def RForce(r): return -(r*np.sin(np.pi*r)+0.5)/S_norm
            def dRForce(r,d_r): return -(r*np.cos(np.pi*r)*(np.pi*d_r) + d_r*np.sin(np.pi*r))/S_norm
        else:
            ValueError("Unsupported force type for numpy disk representation: "+self.force)

        if dfd is None:
            xrot = self.transform(x,x0,HH,yaw,self.dom.Ground)
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
            xrot = self.transform(x,x0,HH,yaw,self.dom.Ground)
            d_xrot = self.transform(x,x0,HH,yaw,self.dom.Ground,dfd=dfd)
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
            tf1  = np.zeros((len(sparse_ids)*dim,1))
            tf2  = np.zeros((len(sparse_ids)*dim,1))
            tf3  = np.zeros((len(sparse_ids)*dim,1))

            ### Fill the output
            tfs = [tf1,tf2,tf3]
            n=[n1,n2,n3]
            for i in range(3):
                tfs[i][0::dim] = d_actuators_x*n[i]
                tfs[i][1::dim] = d_actuators_y*n[i]

        elif dfd == "axial":
            xrot = self.transform(x,x0,HH,yaw,self.dom.Ground)
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
            tf1  = np.zeros((len(sparse_ids)*dim,1))
            tf2  = np.zeros((len(sparse_ids)*dim,1))
            tf3  = np.zeros((len(sparse_ids)*dim,1))

            ### Fill the output
            tfs = [tf1,tf2,tf3]
            n=[n1,n2,n3]
            for i in range(3):
                tfs[i][0::dim] = d_actuators_x*n[i]
                tfs[i][1::dim] = d_actuators_y*n[i]

        elif dfd == "yaw":
            xrot = self.transform(x,x0,HH,yaw,self.dom.Ground)
            d_xrot = self.transform(x,x0,HH,yaw,self.dom.Ground,dfd=dfd)
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
            tf1  = np.zeros((len(sparse_ids)*dim,1))
            tf2  = np.zeros((len(sparse_ids)*dim,1))
            tf3  = np.zeros((len(sparse_ids)*dim,1))

            ### Fill the output
            tfs = [tf1,tf2,tf3]
            n=[n1,n2,n3]
            d_n=[d_n1,d_n2,d_n3]
            for i in range(3):
                tfs[i][0::dim] = (d_actuators_x*n[i]+actuators_x*d_n[i])
                tfs[i][1::dim] = (d_actuators_y*n[i]+actuators_y*d_n[i])

        else:
            raise ValueError("Cannot take the derivative with respect to: "+dfd)

        ### Output the actuator information if needed ### #setting external values this way is dangerous 
        if actuator_disk is not None:
            actuator_array = np.zeros((N*dim,1))
            actuator_array[dim*(sparse_ids)+0] = actuators_x
            actuator_array[dim*(sparse_ids)+1] = actuators_y
            actuator_disk.vector()[:] = actuator_array.T[0]

        return [tf1,tf2,tf3]

    def turbine_force(self,u,inflow_angle,fs):
        
            ### get dof coordinates ###
            x = fs.tf_V0.tabulate_dof_coordinates().T

            ### get only the relevant dofs
            self.sparse_ids = self.construct_sparse_ids(x,1.5)

            ### compute the space kernal and radial force
            self.actuator_disk = Function(fs.tf_V)
            [tf1, tf2, tf3] = self.build_actuator_disk(x,inflow_angle,fs,self.sparse_ids,actuator_disk=self.actuator_disk)

            ### Compose full turbine force
            self.tf = tf1*u[0]**2+tf2*u[1]**2+tf3*u[0]*u[1]

            return self.tf
