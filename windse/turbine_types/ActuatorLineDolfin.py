import numpy as np
import scipy.interpolate as interp
import glob
import os

from windse import windse_parameters
if windse_parameters.dolfin_adjoint:
    from windse.blocks import blockify, ActuatorLineForceBlock

from . import GenericTurbine

from . import Constant, Expression, Function, Point, assemble, dot, dx

from windse.helper_functions import mpi_eval

'''
FIXME After Jeff Meeting 1/13/22

1. figure out how to get rid of u_k in ActuatorLineForceBlock.py (kind of done)
2. take the derivatives with respect to mpi_u_fluid (kind of done)

'''

class ActuatorLine(GenericTurbine):

    def __init__(self, i,x,y,dom,imported_params=None):
        # Define the acceptable column of an wind_farm.csv imput file        
        self.yaml_inputs = ["HH", "RD", "yaw", "yaw", "rpm", "read_turb_data", "blade_segments", "use_local_velocity", "chord_factor", "gauss_factor"]

        # Init turbine
        super(ActuatorLine, self).__init__(i,x,y,dom,imported_params)

        # blockify custom functions so dolfin adjoint can track them
        if self.params.performing_opt_calc:
            block_kwargs = {
                # "construct_sparse_ids":self.construct_sparse_ids,
                "control_types": self.params["optimization"]["control_types"],
                "turb": self
            }
            self.turbine_force = blockify(self.turbine_force,ActuatorLineForceBlock,block_kwargs=block_kwargs)


        # init some flags
        self.first_call_to_alm = True
        self.simTime_prev = None
        self.DEBUGGING = False

        # self.turbine_motion_freq = None#0.1
        # self.turbine_motion_amp = np.radians(20.0)

        self.along_blade_quantities = {}
        self.along_blade_quantities['lift'] = []
        self.along_blade_quantities['drag'] = []
        self.along_blade_quantities['aoa'] = []
        self.along_blade_quantities['axial'] = []
        self.along_blade_quantities['force_x'] = []
        self.along_blade_quantities['force_y'] = []
        self.along_blade_quantities['force_z'] = []

    def get_baseline_chord(self):
        '''
        This function need to return the baseline chord as a numpy array with length num_blade_segments
        '''
        return self.baseline_chord


    def load_parameters(self):
        self.type = self.params["turbines"]["type"]
        self.HH = self.params["turbines"]["HH"]
        self.RD = self.params["turbines"]["RD"]
        self.yaw = self.params["turbines"]["yaw"]
        self.rpm = self.params["turbines"]["rpm"]
        self.read_turb_data = self.params["turbines"]["read_turb_data"]
        self.blade_segments = self.params["turbines"]["blade_segments"]
        self.use_local_velocity = self.params["turbines"]["use_local_velocity"]
        self.chord_factor = self.params["turbines"]["chord_factor"]
        self.gauss_factor = self.params["turbines"]["gauss_factor"]
        self.tip_loss = self.params["turbines"]["tip_loss"]
        self.hub_rad = self.params["turbines"]["hub_rad"]


    def compute_parameters(self):
        # compute turbine radius
        self.radius = 0.5*self.RD
        self.angular_velocity = 2.0*pi*self.rpm/60.0
        self.eps = self.gauss_factor*self.dom.global_hmin


    def create_controls(self, initial_call_from_setup=True):

        if initial_call_from_setup:
            self.controls_list = ["x","y","cl","cd","chord","yaw"] # this is just part of the function as an example of the types of controls 

            self.mx     = Constant(self.x, name="x_{:d}".format(self.index))
            self.my     = Constant(self.y, name="y_{:d}".format(self.index))
            self.myaw   = Constant(self.yaw, name="yaw_{:d}".format(self.index))

        else:
            # Assign the chord, twist, cl, and cd values into a list of constants
            self.mchord = []
            self.mtwist = []
            self.mcl = []
            self.mcd = []

            for k in range(self.num_blade_segments):
                self.mchord.append(Constant(self.chord[k]))
                self.mtwist.append(Constant(self.twist[k]))
                self.mcl.append(Constant(self.cl[k]))
                self.mcd.append(Constant(self.cd[k]))


    def fake_scipy_lift(self, rdim, aoa):
        
        return 1.0


    def fake_scipy_drag(self, rdim, aoa):
        
        return 0.1


    # def mpi_eval(self, u_k, x_0):
    #     a = Constant(u_k(x_0))
        
    #     return a


    def rot_x(self, theta):
        Rx = as_tensor([[1, 0, 0],
                       [0, cos(theta), -sin(theta)],
                       [0, sin(theta), cos(theta)]])
        
        return Rx


    def calculate_relative_fluid_velocity(self, u_k, x_0, n_0, rdim):

        # Evaluate the fluid velocity u_k at the point x_0 
        vel_fluid = mpi_eval(u_k, x_0)
        
        # Remove the component of the fluid velocity oriented along the blade's axis
        vel_fluid -= dot(vel_fluid, n_0[:, 1])*n_0[:, 1]

        # Calculate the relative velocity contribution due to the blade
        # Note, the negative sign is because we seek the velocity of the fluid
        # relative to a stationary blade
        vel_blade = -(rdim+0.01)*self.angular_velocity*n_0[:, 2]

        # The relative velocity is the sum of the fluid and blade velocity
        vel_rel = vel_fluid + vel_blade
        
        # Calculate the magnitude of the relative velocity
        vel_rel_mag = sqrt(vel_rel[0]**2.0 + vel_rel[1]**2.0 + vel_rel[2]**2.0)
        vel_rel_unit = vel_rel/vel_rel_mag

        return vel_rel, vel_rel_mag, vel_rel_unit


    def calcuate_aoa(self, vel_rel, n_0):
            
        def get_angle_between_vectors(a, b, n):
            
            DEBUGGING = False
            
            if DEBUGGING:
                a_x_b = np.dot(np.cross(n, a), b)

                norm_a = np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
                norm_b = np.sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])

                c1 = a_x_b/(norm_a*norm_b)
                c1 = np.clip(c1, -1.0, 1.0)
                aoa_1_np= np.arcsin(c1)

                c2 = np.dot(a, b)/(norm_a*norm_b)
                c2 = np.clip(c2, -1.0, 1.0)
                aoa_2_np = np.arccos(c2)




                a_x_b = dot(cross(n, a), b)

                norm_a = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
                norm_b = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])

                c1 = a_x_b/(norm_a*norm_b)
        #         c1 = np.clip(c1, -1.0, 1.0)
                aoa_1 = asin(c1)

                c2 = dot(a, b)/(norm_a*norm_b)
        #         c2 = np.clip(c2, -1.0, 1.0)
                aoa_2 = acos(c2)

                if aoa_2_np > pi/2.0:
                    if aoa_1_np < 0:
                        aoa_1 = -pi - aoa_1
                    else:
                        aoa_1 = pi - aoa_1
            else:
                aoa_1 = 1.0

            return aoa_1

    #     wind_vec = vel_rel

    #     # Remove the component in the radial direction (along the blade span)
    #     wind_vec -= np.dot(vel_rel, n_0[:, 1])*n_0[:, 1]

    #     vel_fluid -= dot(vel_fluid, n_0[:, 1])*n_0[:, 1]
        # aoa = get_angle_between_vectors(arg1, arg2, arg3)
        # arg1 = in-plane vector pointing opposite rotation (blade sweep direction)
        # arg2 = relative wind vector at node k, including blade rotation effects (wind direction)
        # arg3 = unit vector normal to plane of rotation, in this case, radially along span
        
        # TODO : verify if this is necessary or not
    #     print(vel_rel.values())
    #     vel_rel = vel_rel - dot(vel_rel, n_0[:, 1]) * n_0[:, 1]
    #     print(vel_rel.values())

        aoa = get_angle_between_vectors(-n_0[:, 2], vel_rel, -n_0[:, 1])

        return aoa


    def lookup_lift_coeff(self, rdim, aoa):
        
    #     TODO: Jeff, make this take in rdim and aoa, analogous to looking up terrain 
    #     val = self.dom.Ground(float(x),float(y),dx=dx,dy=dy)
    #     return Constant(val)
        cl = self.fake_scipy_lift(rdim, aoa)

        return Constant(cl)


    def lookup_drag_coeff(self, rdim, aoa):
        
    #     TODO: Jeff, make this take in rdim and aoa, analogous to looking up terrain 
    #     val = self.dom.Ground(float(x),float(y),dx=dx,dy=dy)
    #     return Constant(val)
        cd = self.fake_scipy_drag(rdim, aoa)
        
        return Constant(cd)


    def build_lift_drag_vec(self, vel_rel_unit, n_0):
        # The drag unit vector is oriented opposite of the relative fluid velocity
        drag_unit = -vel_rel_unit
        
        # The lift unit vector is perpendicular to the drag unit vector and the blade axis
        lift_unit = cross(drag_unit, n_0[:, 1])
        
        return lift_unit, drag_unit


    def build_gauss_kernel_at_point(self, u_k, x_0):

        # Get the coordinates of the mesh points as UFL term
        x_mesh = SpatialCoordinate(u_k.function_space().mesh())
        
        # Calculate distance from actuator node, x_0, in 3 directions
        delta_x = x_mesh[0] - x_0[0]
        delta_y = x_mesh[1] - x_0[1]
        delta_z = x_mesh[2] - x_0[2]
        
        # Calculate the distance squared, i.e., (|x_mesh-x|)^2
        r2 = delta_x**2.0 + delta_y**2.0 + delta_z**2.0
        
        # Compute the Gaussian kernel
        gauss_kernel = exp(-r2/(self.eps**2.0))/(self.eps**3.0*np.pi**1.5)

        return gauss_kernel


    def calculate_tip_loss(self, rdim, aoa):
        use_tip_loss = True
    #     if self.tip_loss:
        if use_tip_loss:
            if rdim < 1e-12:
                tip_loss_fac = 1.0

            else:
                loss_exponent = 3.0/2.0*(self.radius-rdim)/(rdim*sin(aoa))
                acos_arg = exp(-loss_exponent)
        #         acos_arg = np.clip(acos_arg, -1.0, 1.0)
                tip_loss_fac = 2.0/pi*acos(acos_arg)

        else:
            tip_loss_fac = 1.0
            
        return tip_loss_fac


    def build_actuator_node(self, u_k, x_0, n_0, rdim):
        """
        Build the force function for a single actuator node.

        This function takes as an input the fluid velocity, location, orientation (n_0),
        and radial dimension of an actuator node along a turbine blade. It returns
        a vector function representing the force as projected onto the computational
        grid. The projection/smoothing step is handled by the spherical Gaussian kernel
        and the magnitude is set by a lift and drag calculation.
        
        Args:
            u_k (FEniCS function):
                A FEniCS/dolfin function representing the current fluid velocity.
            x_0 (``as_vector(1x3 array)``):
                The 3D coordinates specifying the location of the kernel. Note that this
                must reflect the updated location capturing any and all rotation/movement 
                effects.
            n_0 (``as_tensor(3x3 array)``):
                The orientation of the blade at this actuator point.
                n_0[:, 0] = points along the rotor shaft in the +x direction (normal to rotor plane)
                n_0[:, 1] = points along the blade span/axis in the +y direction (along the blade)
                n_0[:, 2] = points tangential to the blade axis (generates a torque about rotor shaft)
            rdim (float):
                The radial dimension in the range [0, blade_radius], with units of meters.

        Returns:
            actuator_force (FEniCS form):
                A vector-valued UFL representation of the turbine force including both
                lift and drag effects for a single actutor point.
        """

        # Calculate the relative fluid velocity, a function of fluid and blade velocity    
        vel_rel, vel_rel_mag, vel_rel_unit = self.calculate_relative_fluid_velocity(u_k, x_0, n_0, rdim)

        # Calculate the angle of attack for this actuator node node
        aoa = self.calcuate_aoa(vel_rel, n_0)
        
        # Lookup the lift and drag coefficients
        cl = self.lookup_lift_coeff(rdim, aoa)
        cd = self.lookup_drag_coeff(rdim, aoa)

        # Calculate the tip loss factor
        tip_loss_fac = self.calculate_tip_loss(rdim, aoa)
                
        # Calculate the magnitude (a scalar quantity) of the lift and drag forces
        lift_mag = tip_loss_fac*(0.5*cl*rho*chord*w*vel_rel_mag**2.0)
        drag_mag = tip_loss_fac*(0.5*cd*rho*chord*w*vel_rel_mag**2.0)
        
        # Build a spherical Gaussian with width eps about point x_0
        gauss_kernel = self.build_gauss_kernel_at_point(u_k, x_0)
        
        # Smear the lift and drag forces onto the grid using the Gaussian kernel
        # The result is a scalar function at each point on the grid
        projected_lift_mag = lift_mag*gauss_kernel
        projected_drag_mag = drag_mag*gauss_kernel
        
        # Compute the orientation of the lift and drag unit vectors
        lift_unit, drag_unit = self.build_lift_drag_vec(vel_rel_unit, n_0)
        
        # The smeared scalar force magnitude is oriented along the lift and drag directions
        lift_force = projected_lift_mag*lift_unit
        drag_force = projected_drag_mag*drag_unit

        # The final force of this actuator is the sum of lift and drag effects
        actuator_force = lift_force + drag_force
            
        return actuator_force


    def turbine_force(self, u, inflow_angle, fs, **kwargs):
        # if dfd is None, alm_output is a dolfin function (tf) [1 x numPts*ndim]
        # otherwise, it returns a numpy array of derivatives [numPts*ndim x numControls]

        try:
            self.simTime = kwargs['simTime']
            self.simTime_prev = kwargs['simTime_prev']
            self.dt = kwargs['dt']
        except:
            raise ValueError('"simTime" and "dt" must be specified for the calculation of ALM force.')


        rho = 1.0

        # for all actuator nodes:
        print('starting actuator node')

        tf = 0
        num_blades = 5
        num_actuator_nodes = 10

        # simTime = Constant(0.0) # simulation time in seconds

        # rpm = 10.0
        # diameter = 130.0
        # radius = 0.5*diameter
        # hub_height = 90


        # hmin = mesh.hmin()/sqrt(3)
        # print(f'hmin = {float(hmin)}')
        # eps = 2.0*hmin

        for blade_id in range(num_blades):
            # This can be built on a blade-by-blade basis
            n_0 = as_tensor([[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]])
            
            theta_0 = 2.0*pi*(blade_id/num_blades)
            theta_0 += self.simTime*self.angular_velocity
            
            rx = self.rot_x(theta_0)
            
            # Why does this need to be transposed to work correctly?
            n_0 = dot(rx, n_0).T
            
            for actuator_id in range(num_actuator_nodes):
                # These quantities may vary on a node-by-node basis
                rdim = self.radius*actuator_id/(num_actuator_nodes-1)
                w = self.radius/(num_actuator_nodes-1)
                
                if rdim < 1e-6 or rdim > self.radius - 1e-6:
                    w = 0.5*w
                    
                chord = Constant(4.0)

                x_0 = as_vector([0.0, rdim, 0.0])
                x_0 = dot(x_0, rx)
                x_0 += as_tensor([0.0, 0.0, hub_height])
                
                x_0_np = np.array(x_0, dtype=float)
                plt.plot(x_0_np[1], x_0_np[2], 'o', color=f'C{blade_id}')
                
                tf += self.build_actuator_node(u_k, x_0, n_0, rdim)

        return tf
        







    def power(self, u, inflow_angle):
        # Create a cylindrical expression aligned with the position of this turbine
        self.cyld_expr = Expression(('sin(yaw)*(x[2]-zs)', '-cos(yaw)*(x[2]-zs)', '(x[1]-ys)*cos(yaw)-(x[0]-xs)*sin(yaw)'),
            degree=1,
            yaw=self.myaw,
            xs=self.mx,
            ys=self.my,
            zs=self.mz)

        self.power_dolfin = assemble(1e-6*dot(-self.tf*self.angular_velocity, self.cyld_expr)*dx)
        self.power_numpy = 1e-6*self.rotor_torque*self.angular_velocity

        # print("in turb.poweer()",self.power_dolfin, self.power_numpy)

        # print('Error between dolfin and numpy: %e' % (float(self.power_dolfin) - self.power_numpy))

        return self.power_dolfin

    def prepare_saved_functions(self, func_list):
        if len(func_list) == 0:
            func_list = [
                [self.tf,"turbine_force"]
            ]
        else:
            func_list[1][0] += self.tf

        return func_list

    def finalize_turbine(self):
        if self.params.rank == 0:
            for key in self.along_blade_quantities.keys():
                data = self.along_blade_quantities[key]
                filename = os.path.join(self.params.folder, f'data/alm/{key}_{self.index}.npy')
                np.save(filename, data)
