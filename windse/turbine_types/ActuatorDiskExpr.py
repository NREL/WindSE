# import windse function
from . import GenericTurbine

# import dolfin functions
from . import Constant, SpatialCoordinate, as_vector, cos, sin, exp, sqrt, dot, assemble, dx

# import windse function
from . import ActuatorDisk

# import other modules
import numpy as np
from scipy.special import gamma

# import from dolfin
from . import Expression

class ActuatorDiskExpr(ActuatorDisk):

    def build_actuator_disk(self,inflow_angle):


        ### Setup some constants ###
        S_norm = (2.0+np.pi)/(2.0*np.pi)
        T_norm = 2.0*gamma(7.0/6.0)
        if self.dom.dim == 3:
            D_norm = np.pi*gamma(4.0/3.0)
        else:
            D_norm = 2.0*gamma(7.0/6.0)

        # The key to getting this to work is to build up the full function
        # using string of smaller parts. 

        ### Create the thrust coeff
        C_tprime = "4*ma/(1-ma)"

        ### Rotate and Shift the Turbine ###
        xrot = "(cos(yaw)*(x[0]-mx) + sin(yaw)*(x[1]-my))"
        yrot = "(-sin(yaw)*(x[0]-mx) + cos(yaw)*(x[1]-my))"
        if self.dom.dim == 3:
            zrot = "(x[2]-mz)"
            A = "(pi*pow(R,2.0))" 
        else:
            A = "(2*R)"
            zrot = "(0.0)"
        xs = [xrot,yrot,zrot]

        ### Create the function that represents the Thickness of the turbine ###
        T = f"(exp(-pow(({xs[0]}/W),6.0)))"

        ### Create the function that represents the Disk of the turbine
        r = f"(sqrt(pow({xs[1]},2.0)+pow({xs[2]},2.0))/R)"
        D = f"(exp(-pow({r},6.0)))"

        ### Create the function that represents the force ###
        if self.force == "constant":
            force = "(1.0)"
        elif self.force == "sine":
            force = f"(({r}*sin(pi*{r})+0.5)/S_norm)"
        else:
            raise NotImplementedError(f"the force type: {self.force}, is not implemented for this turbine type")
        F = f"(-0.5*{A}*{C_tprime}*{force})"

        ## Calculate normalization constant ###
        volNormalization = f"(T_norm*D_norm*W*pow(R,{self.dom.dim-1}))"

        # Build the scalar actuator
        act = f"({F}*{T}*{D}/{volNormalization})"

        # Build the vector actuator
        if self.dom.dim == 3:
            act_list = [f"({act}*cos(yaw))",f"({act}*sin(yaw))",f"(0.0)"]
        else:
            act_list = [f"({act}*cos(yaw))",f"({act}*sin(yaw))"]

        # assemble the full expression from the string fragments:
        actuator_disks = Expression(act_list, 
            degree = 6, 
            domain = self.dom.mesh,
            mx = self.mx, 
            my = self.my, 
            mz = self.mz,
            yaw = self.myaw,
            ma = self.maxial,
            W = self.thickness*1.0,
            R = self.RD/2.0,
            S_norm = S_norm, 
            T_norm = T_norm, 
            D_norm = D_norm, 
            pi = np.pi)

        return actuator_disks

    def turbine_force(self,u,inflow_angle,fs):
        """
        This function creates a turbine force by applying 
        a spacial kernel to each turbine. This kernel is 
        created from the turbines location, yaw, thickness, diameter,
        and force density. Currently, force density is limit to a scaled
        version of 

        .. math::

            r\\sin(r),

        where :math:`r` is the distance from the center of the turbine.

        Args:
            u: the velocity field
            inflow_angle: direction the wind is blowing (0 is west to east and positive rotates clockwise)
            fs: the function space manager object

        Returns:
            tf (dolfin.Function): the turbine force.

        Todo:
            * Setup a way to get the force density from file
        """

        ### compute the space kernal and radial force
        self.inflow_angle.assign(inflow_angle)
        self.actuator_disk = self.build_actuator_disk(self.inflow_angle)

        ### Expand the dot product
        yaw = self.myaw+self.inflow_angle
        tf1 = self.actuator_disk * cos(yaw)**2
        tf2 = self.actuator_disk * sin(yaw)**2
        tf3 = self.actuator_disk * 2.0 * cos(yaw) * sin(yaw)

        ### Compose full turbine force
        self.tf = tf1*u[0]**2+tf2*u[1]**2+tf3*u[0]*u[1]

        return self.tf
