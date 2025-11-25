# import windse function
from typing import Generic
from . import GenericTurbine

# import dolfin functions
from . import Constant, SpatialCoordinate, as_vector, cos, sin, exp, ln, sqrt, dot, assemble, dx

from . import ActuatorDisk

import numpy as np
from scipy.special import gamma

class ActuatorDiskTSPowerCurve(GenericTurbine):

    def __init__(self, i, x, y, dom, imported_params=None):
        # Define the acceptable column of an wind_farm.csv imput file
        self.yaml_inputs = ["HH", "RD", "yaw", "thickness", "CTprime0", "CPprime0", "Vrated", "force"]

        # Init turbine
        super(ActuatorDiskTSPowerCurve, self).__init__(i,x,y,dom,imported_params)

    def load_parameters(self):
        self.HH        = self.params["turbines"]["HH"]
        self.RD        = self.params["turbines"]["RD"]
        self.yaw       = self.params["turbines"]["yaw"]
        self.thickness = self.params["turbines"]["thickness"]
        self.axial     = -1 # needed for generic turbine
        self.CPprime0  = self.params["turbines"]["CPprime0"]
        self.CTprime0  = self.params["turbines"]["CTprime0"]
        self.Prated    = self.params["turbines"]["Prated"]
        self.Trated    = self.params["turbines"]["Trated"]
        self.force     = self.params["turbines"]["force"]

    def create_controls(self):
        self.mx        = Constant(self.x, name="x_{:d}".format(self.index))
        self.my        = Constant(self.y, name="y_{:d}".format(self.index))
        self.myaw      = Constant(self.yaw, name="yaw_{:d}".format(self.index))
        self.mCPprime0 = Constant(self.CPprime0, name="CPprime0_{:d}".format(self.index))
        self.mCTprime0 = Constant(self.CTprime0, name="CTprime0_{:d}".format(self.index))
        self.mPrated   = Constant(self.Prated, name="Prated_{:d}".format(self.index))
        self.mTrated   = Constant(self.Trated, name="Trated_{:d}".format(self.index))

        # The control list is very important for extracting information from the optimzation step.
        self.controls_list = ["x","y","yaw","CPprime0","CTprime0","Prated","Trated"]

    def build_actuator_disk(self):

        ### Alias useful values ###
        x0 = [self.mx,self.my,self.mz]
        yaw = self.myaw+self.inflow_angle
        W = self.thickness*1.0
        R = self.RD/2.0
        CTprime0 = self.mCTprime0
        # CPprime0 = self.mCPprime0
        x=SpatialCoordinate(self.dom.mesh)

        ### Set up some dim dependent values ###
        S_norm = (2.0+np.pi)/(2.0*np.pi)
        T_norm = 2.0*gamma(7.0/6.0)
        if self.dom.dim == 3:
            WTGbase = as_vector((cos(yaw),sin(yaw),0.0))
            A = np.pi*R**2.0
            D_norm = np.pi*gamma(4.0/3.0)
        else:
            WTGbase = as_vector((cos(yaw),sin(yaw)))
            A = 2*R
            D_norm = 2.0*gamma(7.0/6.0)

        ### Rotate and Shift the Turbine ###
        xs = self.yaw_turbine(x,x0,yaw)

        ### Create the function that represents the Thickness of the turbine ###
        T = exp(-pow((xs[0]/W),6.0))

        ### Create the function that represents the Disk of the turbine
        r = sqrt(xs[1]**2.0+xs[2]**2.0)/R
        D = exp(-pow(r,6.0))

        ### Create the function that represents the force ###
        if self.force == "constant":
            force = 1.0
        elif self.force == "sine":
            force = (r*sin(np.pi*r)+0.5)/S_norm
        elif self.force == "chord":
            raise NotImplementedError("not implemented! -cfrontin")
            # self._setup_chord_force()
            # chord = self.mchord[i]
            # force = self.radial_chord_force(r,chord)
        F = -0.5*A*CTprime0*force

        ### Calculate normalization constant ###
        volNormalization = T_norm*D_norm*W*R**(self.dom.dim-1)
        # volNormalization = assemble(T*D*dx)

        # self.fprint(f"Turbine {self.index}:")
        # self.fprint(f"Assembled Volume: {assemble(T*D*dx)}",offset=1)
        # self.fprint(f"Computed Volume:  {volNormalization}",offset=1)
        # self.fprint(f"Assembled Force:  {assemble(F*T*D/volNormalization*dx)}",offset=1)

        # compute disk averaged velocity in yawed case and don't project
        actuator_disks=F*T*D*WTGbase/volNormalization
        # actuator_disks=project(actuator_disks,self.fs.V,solver_type='cg',preconditioner_type="hypre_amg")


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
        self.fs = fs
        self.inflow_angle.assign(inflow_angle)
        self.actuator_disk = self.build_actuator_disk()

        ### Expand the dot product
        yaw = self.myaw+self.inflow_angle
        self.tf1 = self.actuator_disk * cos(yaw)**2
        self.tf2 = self.actuator_disk * sin(yaw)**2
        self.tf3 = self.actuator_disk * 2.0 * cos(yaw) * sin(yaw)

        # compute power curve adjustment as a function of velocity
        CPprime0 = self.mCPprime0
        CTprime0 = self.mCTprime0
        Prated = self.mPrated
        Trated = self.mTrated
        Vrated = Prated*CTprime0/(Trated*CPprime0)
        A = np.pi/4.0*self.RD**2.0
        Vps = sqrt(2*Prated/(A*Vrated*CPprime0))

        # precompute key values
        beta_smooth = 16.0
        vel_magnitude = sqrt(u[0]**2 + u[1]**2 + u[2]**2)

        f0 = (0.5*A*vel_magnitude**3)/1e6
        f1 = (0.5*A*vel_magnitude**3)*CTprime0/1e6
        f2 = (0.5*A*vel_magnitude**3)*CTprime0*Vps**2/(vel_magnitude**2)/1e6
        f3 = (0.5*A*vel_magnitude**3)*CTprime0*Vps**2*Vrated/(vel_magnitude**3)/1e6

        # smooth twice: one in power space and one in coefficient space

        # # first attempt: boltzmann operator
        # # first blend: smoothmin Region II.5 and Region III power
        # blend1 = 1.0*(
        #     f2*exp(-beta_smooth*f2) + f3*exp(-beta_smooth*f3)
        # )/(
        #     exp(-beta_smooth*f2) + exp(-beta_smooth*f3)
        # )
        #
        # # second blend: smoothmin Region II and first blend coefficient
        # CTp_factor = (
        #     f1/f0*exp(-beta_smooth*f1/f0) + blend1/f0*exp(-beta_smooth*blend1/f0)
        # )/(
        #     exp(-beta_smooth*f1/f0) + exp(-beta_smooth*blend1/f0)
        # )

        # second attempt: mellowmax operator
        # first blend: smoothmin Region II.5 and Region III power
        blend1 = -1.0/beta_smooth*ln(exp(-beta_smooth*f2) + exp(-beta_smooth*f3))

        # second blend: smoothmin Region II and first blend coefficient
        CTp_factor = -1.0/beta_smooth*ln(
            exp(-beta_smooth*f1/f0) + exp(-beta_smooth*blend1/f0)
        )/CTprime0

        ### Compose full turbine force
        self.tf = CTp_factor*(self.tf1*u[0]**2+self.tf2*u[1]**2+self.tf3*u[0]*u[1])

        return self.tf

    def force_gradient(self):
        pass

    def power(self, u, inflow_angle):
        # adjust for turbine inefficiency
        return self.mCPprime0/self.mCTprime0*dot(-self.tf,u)/1.0e6

    def power_gradient(self):
        pass

    def prepare_saved_functions(self, func_list):
        if len(func_list) == 0:
            func_list = [
                [self.actuator_disk,"actuator_disk"],
                [self.tf,"turbine_force"]
            ]
        else:
            func_list[0][0] += self.actuator_disk
            func_list[1][0] += self.tf

        return func_list
