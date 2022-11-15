import numpy as np
import scipy.interpolate as interp
import glob
import os

from windse import windse_parameters
if windse_parameters.dolfin_adjoint:
    from windse.blocks import blockify, InterpBlock
    from pyadjoint import AdjFloat

from . import GenericTurbine

from . import (Constant, Expression, Function, Point, assemble, dot, dx,
pi, cos, acos, asin, sin, sqrt, exp, cross, as_tensor, as_vector, SpatialCoordinate,
project, inner, TrialFunction, TestFunction, PETScKrylovSolver, MPI)
from windse.helper_functions import mpi_eval, ufl_eval, test_dolfin_adjoint

'''
FIXME After Jeff Meeting 1/13/22

1. figure out how to get rid of u_k in ActuatorLineForceBlock.py (kind of done)
2. take the derivatives with respect to mpi_u_fluid (kind of done)

'''

class ActuatorLineDolfin(GenericTurbine):

    def __init__(self, i,x,y,dom,imported_params=None):
        # Define the acceptable column of an wind_farm.csv imput file        
        self.yaml_inputs = ["HH", "RD", "yaw", "yaw", "rpm", "read_turb_data", "blade_segments", "use_local_velocity", "chord_factor", "gauss_factor"]

        # Init turbine
        super(ActuatorLineDolfin, self).__init__(i,x,y,dom,imported_params)

        # blockify custom functions so dolfin adjoint can track them
        if self.params.performing_opt_calc:
            block_kwargs = {
                "func": self.lookup_lift_coeff
            }
            self.lookup_lift_coeff = blockify(self.lookup_lift_coeff,InterpBlock,block_kwargs=block_kwargs)
            block_kwargs = {
                "func": self.lookup_drag_coeff
            }
            self.lookup_drag_coeff = blockify(self.lookup_drag_coeff,InterpBlock,block_kwargs=block_kwargs)


        # init some flags
        self.first_call_to_alm = True
        self.simTime_prev = None
        self.DEBUGGING = False

        # self.turbine_motion_freq = None#0.1
        # self.turbine_motion_amp = np.radians(20.0)

    def get_baseline_chord(self):
        '''
        This function need to return the baseline chord as a numpy array with length num_actuator_nodes
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
        self.chord_perturb = float(self.params["turbines"]["chord_perturb"])
        self.chord_perturb_id = int(self.params["turbines"]["chord_perturb_id"])
        self.chord_override = self.params["turbines"]["chord_override"]
        self.motion_file = self.params["turbines"]["motion_file"]
        self.use_gauss_vel_probe = self.params["turbines"]["use_gauss_vel_probe"]


    def compute_parameters(self):
        pass

    def create_controls(self, initial_call_from_setup=True):

        if initial_call_from_setup:
            self.controls_list = ["x","y","chord","yaw","twist"] # this is just part of the function as an example of the types of controls 

            self.mx     = Constant(self.x, name="x_{:d}".format(self.index))
            self.my     = Constant(self.y, name="y_{:d}".format(self.index))
            self.myaw   = Constant(self.yaw, name="yaw_{:d}".format(self.index))

        else:
            # Assign the chord, twist, cl, and cd values into a list of constants
            self.mchord = []
            self.mtwist = []
            self.mcl = []
            self.mcd = []

            for i in range(self.num_blades):
                self.mcl.append([])
                self.mcd.append([])
                for j in range(self.num_actuator_nodes):
                    self.mcl[i].append(Constant(self.cl[j],name=f"cl_{i,j}"))
                    self.mcd[i].append(Constant(self.cd[j],name=f"cd_{i,j}"))

            for j in range(self.num_actuator_nodes):
                self.mchord.append(Constant(self.chord[j],name=f"chord_{j}"))
                self.mtwist.append(Constant(self.twist[j],name=f"twist_{j}"))

    def init_platform_motion(self):
        # 0: Time (s)
        # 1: PtfmSurge (m)
        # 2: PtfmSway (m)
        # 3: PtfmHeave (m)
        # 4: PtfmRoll (deg)
        # 5: PtfmPitch (deg)
        # 6: PtfmYaw (deg)

        path_to_motion_file = os.path.join(os.path.dirname(self.read_turb_data), self.motion_file)
        motion_data = np.genfromtxt(path_to_motion_file, skip_header=1)

        time_0 = motion_data[:, 0]
        pitch_0 = motion_data[:, 5]

        self.motion_interp = interp.interp1d(time_0, pitch_0, kind='linear')

    def calc_platform_motion(self):

        if hasattr(self, 'motion_interp'):
            platform_theta = self.motion_interp(float(self.simTime))
            platform_theta = np.radians(platform_theta)
            # print(platform_theta)

        else:
            platform_theta = 0.0

        return platform_theta


        '''
        pos_before_floating = self.rotate_points(pos, np.radians(wave_interp(theta_prev_tt)), [0, 1, 0])
        pos, unit_vec, vel = self.rotate_points([pos, unit_vec, vel], np.radians(wave_interp(theta_tt)), [0, 1, 0])
        pos_prev = self.rotate_points(pos_prev, np.radians(wave_interp(theta_prev_tt)), [0, 1, 0])

        # Finite different calc of the velocity incurred by this shift in the floating position
        # Negative because we seek the velocity of fluid relative to a stationary blade
        if self.simTime_prev is None:
            vel_due_to_floating = 0.0
        else:
            vel_due_to_floating = -(pos - pos_before_floating)/(0.5*(self.simTime + self.dt - self.simTime_prev))

        test_experimental_finite_diff = True

        if test_experimental_finite_diff:
            vel += vel_due_to_floating
        '''

    def init_alm_calc_terms(self):
        # compute turbine radius
        self.rho = 1.0
        self.num_blades = 3
        self.radius = 0.5*self.RD
        self.angular_velocity = 2.0*pi*self.rpm/60.0


        self.simTime = Constant(0.0,name="simTime")
        self.simTime_prev = Constant(0.0,name="simTime_prev")
        self.dt = Constant(0.0,name="dt")

        self.eps = self.gauss_factor*self.dom.global_hmin

        if self.blade_segments == "computed":
            self.num_actuator_nodes = int(2.0*self.radius/self.eps)
        else:
            self.num_actuator_nodes = self.blade_segments


        place_nodes_at_extremes = False

        if place_nodes_at_extremes:
            w = self.radius/(self.num_actuator_nodes-1)
            self.w = []
            self.rdim = []
            for i in range(self.num_actuator_nodes):
                if i == 0 or i == self.num_actuator_nodes - 1:
                    self.w.append(Constant(w/2.0, name=f"w_{i}"))
                else:
                    self.w.append(Constant(w, name=f"w_{i}"))
                self.rdim.append(Constant(self.radius*i/(self.num_actuator_nodes-1),name=f"rdim_{i}"))

        else:
            w = self.radius/self.num_actuator_nodes
            self.w = []
            self.rdim = []
            for i in range(self.num_actuator_nodes):
                self.w.append(Constant(w, name=f"w_{i}"))
                self.rdim.append(Constant(self.radius*(i+0.5)/self.num_actuator_nodes,name=f"rdim_{i}"))


        # Initialize lists
        self.gauss_kernel = []
        self.aoa_forms = []
        self.aoa_values = []
        self.vel_fluid = []
        self.x_0 = []
        self.x_0_prev = []
        self.x_0_pre_motion = []
        
        self.lift_mag = []
        self.drag_mag = []
        self.axial = []
        self.actuator_force_components = []

        for i in range(self.num_blades):
            self.gauss_kernel.append([])
            self.aoa_forms.append([])
            self.aoa_values.append([])
            self.vel_fluid.append([])
            self.x_0.append([])
            self.x_0_prev.append([])
            self.x_0_pre_motion.append([])

            self.lift_mag.append([])
            self.drag_mag.append([])
            self.axial.append([])
            self.actuator_force_components.append([])

            for j in range(self.num_actuator_nodes):
                self.gauss_kernel[i].append(0.0)
                self.aoa_forms[i].append(0.0)
                self.aoa_values[i].append(Constant(0.0,name=f"aoa_{(i,j)}"))
                self.vel_fluid[i].append(as_vector((Constant(0.0,name=f"ux_{(i,j)}"),Constant(0.0,name=f"uy_{(i,j)}"),Constant(0.0,name=f"uz_{(i,j)}"))))
                # self.x_0_prev[i].append(as_vector((Constant(0.0),Constant(0.0),Constant(0.0))))
                self.x_0[i].append(as_vector([0.0, 0.0, 0.0]))
                self.x_0_prev[i].append(as_vector([0.0, 0.0, 0.0]))
                self.x_0_pre_motion[i].append(as_vector([0.0, 0.0, 0.0]))

                self.lift_mag[i].append(Constant(0.0,name=f"lift_mag_{(i,j)}"))
                self.drag_mag[i].append(Constant(0.0,name=f"drag_mag_{(i,j)}"))
                self.axial[i].append(Constant(0.0,name=f"axial_{(i,j)}"))
                self.actuator_force_components[i].append([])

    def init_blade_properties(self):
        if self.read_turb_data:

            self.fprint('Setting chord, lift, and drag from file \'%s\'' % (self.read_turb_data))

            actual_turbine_data = np.genfromtxt(self.read_turb_data, delimiter = ',', skip_header = 1)

            actual_x = actual_turbine_data[:, 0]

            actual_chord = self.chord_factor*actual_turbine_data[:, 1]

            # Baseline twist is expressed in degrees, convert to radians
            actual_twist = actual_turbine_data[:, 2]/180.0*np.pi

            actual_cl = actual_turbine_data[:, 3]
            actual_cd = actual_turbine_data[:, 4]

            # Create interpolators for chord, lift, and drag
            inter_chord = interp.interp1d(actual_x, actual_chord)
            interp_twist = interp.interp1d(actual_x, actual_twist)
            interp_cl = interp.interp1d(actual_x, actual_cl)
            interp_cd = interp.interp1d(actual_x, actual_cd)

            # Construct the points at which to generate interpolated values
            interp_points = np.linspace(0.0, 1.0, self.num_actuator_nodes)

            # Generate the interpolated values
            chord = inter_chord(interp_points)
            twist = interp_twist(interp_points)

            # Note that this is not the lookup table for cl and cd based on AoA and radial
            # location along the blade, rather, it is the initialization based on a
            # 1D interpolation where AoA is assumed to be a few degrees set back from stall.
            # If using a lookup table which include AoA, these values will be overwritten 
            # during the course of the solve.
            cl = interp_cl(interp_points)
            cd = interp_cd(interp_points)

        else:
            # If not reading from a file, prescribe dummy values
            chord = self.radius/20.0*np.ones(self.num_actuator_nodes)
            twist = np.zeros(self.num_actuator_nodes)
            cl = np.ones(self.num_actuator_nodes)
            cd = 0.1*np.ones(self.num_actuator_nodes)

        if self.chord_override is not None:
            chord_new = np.genfromtxt(self.chord_override, delimiter=',')
            assert np.size(chord_new) == np.size(chord)
            chord = np.copy(chord_new)

        if abs(self.chord_perturb) > 0:
            self.fprint(f'Perturbing chord index {self.chord_perturb_id} by {self.chord_perturb}')
            self.fprint(f'Old Value: {chord[self.chord_perturb_id]}')
            self.fprint(f'New Value: {chord[self.chord_perturb_id]+self.chord_perturb}')
            chord[self.chord_perturb_id] += self.chord_perturb

        if self.params.rank == 0:
            chord_filename = os.path.join(self.params.folder, f'data/alm/chord.csv')
            np.savetxt(chord_filename, chord, delimiter=',')

        # Store the chord, twist, cl, and cd values
        self.chord = chord
        self.twist = twist
        self.cl = cl
        self.cd = cd

        # Store a copy of the original values for future use
        self.baseline_chord = chord
        self.baseline_twist = twist
        self.baseline_cl = cl
        self.baseline_cd = cd

        self.max_chord = 1.5*np.amax(chord)


    def init_lift_drag_lookup_tables(self):

        if self.read_turb_data:

            airfoil_data_path = os.path.dirname(self.read_turb_data)+'/airfoil_polars'

            # Determine the number of files in airfoil_data_path
            num_stations = len(glob.glob('%s/*.txt' % (airfoil_data_path)))

            station_radii = np.linspace(0.0, self.radius, num_stations)

            test_min_d_theta = np.zeros(num_stations)
            test_min_angle = np.zeros(num_stations)
            test_max_angle = np.zeros(num_stations)

            for station_id in range(num_stations):
                data = np.genfromtxt('%s/af_station_%d.txt' % (airfoil_data_path, station_id), skip_header=1, delimiter=' ')
                angles = data[:, 0]

                test_min_d_theta[station_id] = np.amin(angles[1:] - angles[:-1])
                test_min_angle[station_id] = np.amin(angles)
                test_max_angle[station_id] = np.amax(angles)

            min_d_theta = np.amin(test_min_d_theta)
            min_angle = np.amin(test_min_angle)
            max_angle = np.amax(test_max_angle)

            if min_d_theta < np.radians(0.5):
                min_d_theta = np.radians(0.5)

            ni = int((max_angle - min_angle)/(0.5*min_d_theta))
            angles_i = np.linspace(min_angle, max_angle, ni)

            for station_id in range(num_stations):
                # print('Reading Airfoil Data #%d' % (station_id))
                data = np.genfromtxt('%s/af_station_%d.txt' % (airfoil_data_path, station_id), skip_header=1, delimiter=' ')

                angles_0 = data[:, 0]
                c_lift_0 = data[:, 1]
                c_drag_0 = data[:, 2]

                s_i = station_radii[station_id]*np.ones(np.size(angles_i))

                c_lift_interp = interp.interp1d(angles_0, c_lift_0, kind='linear')
                c_drag_interp = interp.interp1d(angles_0, c_drag_0, kind='linear')

                c_lift_i = c_lift_interp(angles_i)
                c_drag_i = c_drag_interp(angles_i)

                if station_id == 0:
                    station = s_i
                    angles = angles_i
                    c_lift = c_lift_i
                    c_drag = c_drag_i
                else:
                    station = np.vstack((station, s_i))
                    angles = np.vstack((angles, angles_i))
                    c_lift = np.vstack((c_lift, c_lift_i))
                    c_drag = np.vstack((c_drag, c_drag_i))

            # Create interpolation functions for lift and drag based on angle of attack and location along blade
            if 'linear_interp' in self.params['general']['name']:
                print(f'Using linear interpolator: LinearNDInterpolator')
                nodes = np.vstack((station.flatten(), angles.flatten())).T

                # Create interpolation functions for lift and drag based on angle of attack and location along blade
                self.lookup_cl = interp.LinearNDInterpolator(nodes, c_lift.flatten())
                self.lookup_cd = interp.LinearNDInterpolator(nodes, c_drag.flatten())
                self.lookup_takes_derivatives = False

            else:
                self.lookup_cl = interp.RectBivariateSpline(station_radii, angles_i, c_lift)
                self.lookup_cd = interp.RectBivariateSpline(station_radii, angles_i, c_drag)
                self.lookup_takes_derivatives = True

        else:
            self.lookup_cl = None # TODO: We probably need to handle this case
            self.lookup_cd = None


    def lookup_lift_coeff(self, rdim, aoa, dx=0, dy=0):

        if self.lookup_takes_derivatives:
            cl = self.lookup_cl(rdim, aoa, dx=dx, dy=dy)[0][0]
            # print(f"rank: {MPI.comm_world.Get_rank()}, {float(rdim)}, {float(aoa)}, cl: {cl}")
        else:
            cl = self.lookup_cl(rdim, aoa)


        return cl


    def lookup_drag_coeff(self, rdim, aoa, dx=0, dy=0):
        
        if self.lookup_takes_derivatives:
            cd = self.lookup_cd(rdim, aoa, dx=dx, dy=dy)[0][0]
            # print(f"rank: {MPI.comm_world.Get_rank()}, {float(rdim)}, {float(aoa)}, cd: {cd}")
        else:
            cd = self.lookup_cd(rdim, aoa)

        return cd

    # TODO: use the one from GenericTurbine
    def rot_x(self, theta):
        # Rotation about x-axis (e.g., due to normal rotor rotation)
        Rx = as_tensor([[1,          0,           0],
                        [0, cos(theta), -sin(theta)],
                        [0, sin(theta),  cos(theta)]])

        return Rx

    # TODO: use the one from GenericTurbine
    def rot_y(self, theta):
        # Rotation about y-axis (e.g., due to pitching platform)
        Ry = as_tensor([[ cos(theta), 0, sin(theta)],
                        [          0, 1,          0],
                        [-sin(theta), 0, cos(theta)]])

        return Ry

    # TODO: use the one from GenericTurbine
    def rot_z(self, theta):
        # Rotation about z-axis (e.g., due to yawing)
        Rz = as_tensor([[cos(theta), -sin(theta), 0],
                        [sin(theta),  cos(theta), 0],
                        [         0,           0, 1]])

        return Rz


    def calculate_relative_fluid_velocity(self, u_k, x_0, x_0_pre_motion, n_0, blade_id, actuator_id):
        rdim = self.rdim[actuator_id]

        # vel_fluid_temp = mpi_eval(u_k, x_0)
        # # vel_fluid_temp = [AdjFloat(8.0),AdjFloat(0.0),AdjFloat(0.0)]#mpi_eval(u_k, x_0)
        # for i in range(self.dom.dim):
        #     self.vel_fluid[blade_id][actuator_id][i].assign(vel_fluid_temp[i])
        vel_fluid = self.vel_fluid[blade_id][actuator_id]

        self.axial[blade_id][actuator_id] = dot(vel_fluid, n_0[:, 0])

        # Remove the component of the fluid velocity oriented along the blade's axis
        vel_fluid -= dot(vel_fluid, n_0[:, 1])*n_0[:, 1]

        # Calculate the relative velocity contribution due to the blade
        # Note, the negative sign is because we seek the velocity of the fluid
        # relative to a stationary blade
        vel_blade = -(rdim+0.01)*self.angular_velocity*n_0[:, 2]

        # Calculate the velocity due to any platform motion
        # Note: this finite difference is not as accurate as the analytic method
        # using angular velocity above, so the former should be used if no
        # complicated platform motion is expected.
        # As before, take the negative sign since we seek the velocity of
        # the fluid as relative to a stationary blade.
        if self.motion_file is not None:
            vel_blade_2 = -(x_0 - x_0_pre_motion)/(self.simTime - self.simTime_prev + 1e-6)
            # x_0_np = np.array(x_0, dtype=float)
            # x_0_pre_motion_np = np.array(x_0_pre_motion, dtype=float)
            # vel_blade_2_np = -(x_0_np - x_0_pre_motion_np)/(self.simTime - self.simTime_prev)
            # print('vel_blade = ', np.array(vel_blade_2_np, dtype=float))
        else:
            vel_blade_2 = as_vector([0, 0, 0])

        # The relative velocity is the sum of the fluid and blade velocity
        vel_rel = vel_fluid + vel_blade + vel_blade_2
        
        # Calculate the magnitude of the relative velocity
        vel_rel_mag = sqrt(vel_rel[0]**2.0 + vel_rel[1]**2.0 + vel_rel[2]**2.0)
        vel_rel_unit = vel_rel/vel_rel_mag



        return vel_rel, vel_rel_mag, vel_rel_unit


    def calcuate_aoa(self, vel_rel, n_0, twist):
            
        def get_angle_between_vectors(a, b, n):
            a_x_b = dot(cross(n, a), b)

            norm_a = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
            norm_b = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])

            c1 = a_x_b/(norm_a*norm_b)
            # c1 = np.clip(c1, -1.0, 1.0)
            aoa_1 = asin(c1)

            c2 = dot(a, b)/(norm_a*norm_b)
            # c2 = np.clip(c2, -1.0, 1.0)
            # print(float(c2))
            aoa_2 = acos(c2)

            if float(aoa_2) > pi/2.0:
                if float(aoa_1) < 0:
                    aoa_1 = -pi - aoa_1
                else:
                    aoa_1 = pi - aoa_1

            return aoa_1

        # wind_vec = vel_rel

        # # Remove the component in the radial direction (along the blade span)
        # wind_vec -= np.dot(vel_rel, n_0[:, 1])*n_0[:, 1]

        # # aoa = get_angle_between_vectors(arg1, arg2, arg3)
        # # arg1 = in-plane vector pointing opposite rotation (blade sweep direction)
        # # arg2 = relative wind vector at node k, including blade rotation effects (wind direction)
        # # arg3 = unit vector normal to plane of rotation, in this case, radially along span
        
        # # TODO : verify if this is necessary or not
        # print(vel_rel.values())
        # vel_rel = vel_rel - dot(vel_rel, n_0[:, 1]) * n_0[:, 1]
        # print(vel_rel.values())

        aoa = get_angle_between_vectors(-n_0[:, 2], vel_rel, -n_0[:, 1])

        aoa -= twist

        return aoa

    def build_lift_drag_vec(self, vel_rel_unit, n_0):
        # The drag unit vector is oriented opposite of the relative fluid velocity
        drag_unit = -vel_rel_unit
        
        # The lift unit vector is perpendicular to the drag unit vector and the blade axis
        lift_unit = cross(drag_unit, n_0[:, 1])
        
        return lift_unit, drag_unit


    def build_gauss_kernel_at_point(self, x_0):

        # Get the coordinates of the mesh points as UFL term
        x_mesh = SpatialCoordinate(self.fs.V.mesh())
        
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
        if self.tip_loss:
            if float(rdim) < 1e-12:
                tip_loss_fac = 1.0

            else:
                loss_exponent = 3.0/2.0*(self.radius-rdim)/(rdim*sin(aoa))
                acos_arg = exp(-loss_exponent)
        #         acos_arg = np.clip(acos_arg, -1.0, 1.0)
                tip_loss_fac = 2.0/pi*acos(acos_arg)

        else:
            tip_loss_fac = 1.0
            
        return tip_loss_fac


    def vel_fluid_eval(self,u,g):
        out = []
        for i in range(self.dom.dim):
            out.append(assemble(u[i]*g*dx)/assemble(g*dx))
        return out


    def build_actuator_node(self, u_k, x_0, x_0_pre_motion, n_0, blade_id, actuator_id):
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
        # Get blade properties
        rdim = self.rdim[actuator_id]
        w = self.w[actuator_id]
        twist = self.mtwist[actuator_id]
        chord = self.mchord[actuator_id]

        # Build a spherical Gaussian with width eps about point x_0
        self.gauss_kernel[blade_id][actuator_id] = self.build_gauss_kernel_at_point(x_0)


        # approximate the fluid velocity u_k at the point x_0 
        if self.use_gauss_vel_probe:
            vel_fluid_temp = self.vel_fluid_eval(u_k, self.gauss_kernel[blade_id][actuator_id])
        else:
            vel_fluid_temp = mpi_eval(u_k, x_0)

        for i in range(self.dom.dim):
            self.vel_fluid[blade_id][actuator_id][i].assign(vel_fluid_temp[i])
        vel_fluid = self.vel_fluid[blade_id][actuator_id]

        # Calculate the relative fluid velocity, a function of fluid and blade velocity    
        vel_rel, vel_rel_mag, vel_rel_unit = self.calculate_relative_fluid_velocity(u_k, x_0, x_0_pre_motion, n_0, blade_id, actuator_id)

        # Calculate the angle of attack for this actuator node node
        self.aoa_forms[blade_id][actuator_id] = self.calcuate_aoa(vel_rel, n_0, twist)
        self.aoa_values[blade_id][actuator_id].assign(ufl_eval(self.aoa_forms[blade_id][actuator_id],print_statement=f"aoa: {(blade_id,actuator_id)}, time: {float(self.simTime):1.4f}, rank: {self.params.rank}"))
        aoa_f = self.aoa_values[blade_id][actuator_id]
        
        # Lookup the lift and drag coefficients
        self.mcl[blade_id][actuator_id].assign(self.lookup_lift_coeff(rdim, aoa_f))
        self.mcd[blade_id][actuator_id].assign(self.lookup_drag_coeff(rdim, aoa_f))
        cl = self.mcl[blade_id][actuator_id]
        cd = self.mcd[blade_id][actuator_id]

        # Calculate the tip loss factor
        tip_loss_fac = self.calculate_tip_loss(rdim, self.aoa_forms[blade_id][actuator_id])
                
        # Calculate the magnitude (a scalar quantity) of the lift and drag forces
        lift_mag = tip_loss_fac*(0.5*cl*self.rho*chord*w*vel_rel_mag*vel_rel_mag)
        drag_mag = tip_loss_fac*(0.5*cd*self.rho*chord*w*vel_rel_mag*vel_rel_mag)

        self.lift_mag[blade_id][actuator_id] = lift_mag
        self.drag_mag[blade_id][actuator_id] = drag_mag

        # if actuator_id == 2:
        #     control_list = []
        #     control_list.append(twist)
        #     control_list.append(chord)
        #     test_dolfin_adjoint(control_list,lift_mag)
        #     print(ufl_eval(tip_loss_fac*(0.5*cl*self.rho*w*vel_rel_mag**2.0)))
        #     exit()

        
        # Smear the lift and drag forces onto the grid using the Gaussian kernel
        # The result is a scalar function at each point on the grid
        projected_lift_mag = lift_mag*self.gauss_kernel[blade_id][actuator_id]
        projected_drag_mag = drag_mag*self.gauss_kernel[blade_id][actuator_id]

        # Compute the orientation of the lift and drag unit vectors
        lift_unit, drag_unit = self.build_lift_drag_vec(vel_rel_unit, n_0)
        
        # The smeared scalar force magnitude is oriented along the lift and drag directions
        lift_force = projected_lift_mag*lift_unit
        drag_force = projected_drag_mag*drag_unit

        # The final force of this actuator is the sum of lift and drag effects
        actuator_force = lift_force + drag_force
        # actuator_force = project(actuator_force,self.fs.V, solver_type='cg', preconditioner_type='jacobi')

        actuator_force_components = []

        for i in range(self.dom.dim):
            f_i = dot(actuator_force, n_0[:, i])
            actuator_force_components.append(f_i)

        self.actuator_force_components[blade_id][actuator_id] = actuator_force_components

        return actuator_force


    def turbine_force(self, u, inflow_angle, fs, **kwargs):
        # if dfd is None, alm_output is a dolfin function (tf) [1 x numPts*ndim]
        # otherwise, it returns a numpy array of derivatives [numPts*ndim x numControls]

        if self.first_call_to_alm:
            self.fs = fs
            self.init_alm_calc_terms()
            self.init_blade_properties()
            self.init_lift_drag_lookup_tables()
            self.create_controls(initial_call_from_setup=False)

            if self.motion_file is not None:
                self.init_platform_motion()

            self.fprint(f'Gaussian Width: {self.eps}')
            self.fprint(f'Minimum Mesh Spacing: {self.dom.global_hmin}')
            self.fprint(f'Number of Actuators per Blade: {self.num_actuator_nodes}')


            debugging_chunked_form = True

            if self.num_actuator_nodes > 20 or debugging_chunked_form == True:
                # To avoid running out of memory for what will be a very
                # large UFL form, store the form in multiple chunks associated
                # with each of the turbine blades in a [1 x num_blades] list. 
                # e.g., tf[0] = the UFL form for blade number 1, and so on.
                tf = [0 for k in range(self.num_blades)]
            else:
                tf = 0

            self.along_blade_quantities = {}
            self.along_blade_quantities['lift'] = []
            self.along_blade_quantities['drag'] = []
            self.along_blade_quantities['aoa'] = []
            self.along_blade_quantities['axial'] = []

            self.along_blade_quantities['force_x'] = []
            self.along_blade_quantities['force_y'] = []
            self.along_blade_quantities['force_z'] = []

            along_blade_lift = []
            along_blade_drag = []
            along_blade_aoa = []
            along_blade_axial = []

            along_blade_force_x = []
            along_blade_force_y = []
            along_blade_force_z = []

            platform_theta = self.calc_platform_motion()
            self.platform_theta = Constant(platform_theta, name="platform_theta")
            self.platform_theta_prev = Constant(platform_theta, name="platform_theta_prev")
            
            for blade_id in range(self.num_blades):
                # This can be built on a blade-by-blade basis
                n_0_base = as_tensor([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]])
                
                theta_0 = 2.0*pi*(blade_id/self.num_blades)
                theta = theta_0 + self.simTime*self.angular_velocity
                

                rx = self.rot_x(-theta)


                ry = self.rot_y(-self.platform_theta)
                
                # Why does this need to be transposed to work correctly?
                # n_0 = dot(rx, n_0_base).T # WORKS FOR NO PLATFORM MOTION
                n_0 = dot(ry, dot(rx, n_0_base)).T

                theta_prev = theta_0 + self.simTime_prev*self.angular_velocity
                
                rx_prev = self.rot_x(-theta_prev)
                ry_prev = self.rot_y(-self.platform_theta_prev)

                # Why does this need to be transposed to work correctly?
                # n_0_prev = dot(rx_prev, n_0_base).T # WORKS FOR NO PLATFORM MOTION
                n_0_prev = dot(ry_prev, dot(rx_prev, n_0_base)).T

                for actuator_id in range(self.num_actuator_nodes):

                    x_0_base = as_tensor([0.0, self.rdim[actuator_id], 0.0])

                    self.x_0[blade_id][actuator_id] = dot(x_0_base, rx)
                    self.x_0[blade_id][actuator_id] += as_tensor([0.0, 0.0, self.mz]) # TODO: what about mx and my do those matter
                    self.x_0[blade_id][actuator_id] = dot(self.x_0[blade_id][actuator_id], ry)

                    self.x_0_pre_motion[blade_id][actuator_id] = dot(x_0_base, rx)
                    self.x_0_pre_motion[blade_id][actuator_id] += as_tensor([0.0, 0.0, self.mz]) # TODO: what about mx and my do those matter
                    self.x_0_pre_motion[blade_id][actuator_id] = dot(self.x_0_pre_motion[blade_id][actuator_id], ry_prev)

                    self.x_0_prev[blade_id][actuator_id] = dot(x_0_base, rx_prev)
                    self.x_0_prev[blade_id][actuator_id] += as_tensor([0.0, 0.0, self.mz])
                    self.x_0_prev[blade_id][actuator_id] = dot(self.x_0_prev[blade_id][actuator_id], ry_prev)

                    # tf += self.build_actuator_node(u, self.x_0[blade_id][actuator_id], self.x_0_pre_motion[blade_id][actuator_id], n_0, blade_id, actuator_id)
                    if isinstance(tf, list):
                        tf[blade_id] += self.build_actuator_node(u, self.x_0[blade_id][actuator_id], self.x_0_pre_motion[blade_id][actuator_id], n_0, blade_id, actuator_id)
                    else:
                        tf += self.build_actuator_node(u, self.x_0[blade_id][actuator_id], self.x_0_pre_motion[blade_id][actuator_id], n_0, blade_id, actuator_id)

                    along_blade_lift.append(float(self.lift_mag[blade_id][actuator_id]))
                    along_blade_drag.append(float(self.drag_mag[blade_id][actuator_id]))
                    along_blade_aoa.append(float(self.aoa_values[blade_id][actuator_id]))
                    along_blade_axial.append(float(self.axial[blade_id][actuator_id]))

                    # TODO : wrap the appended quantity (e.g., fx, fy, fz) as a float and store the evaluated value
                    # TODO: wrap these as assemble and only save when necessary, i.e., if saving_blade_loads ...
                    along_blade_force_x.append(self.actuator_force_components[blade_id][actuator_id][0])
                    along_blade_force_y.append(self.actuator_force_components[blade_id][actuator_id][1])
                    along_blade_force_z.append(self.actuator_force_components[blade_id][actuator_id][2])

            self.along_blade_quantities['lift'].append(along_blade_lift)
            self.along_blade_quantities['drag'].append(along_blade_drag)
            self.along_blade_quantities['aoa'].append(along_blade_aoa)
            self.along_blade_quantities['axial'].append(along_blade_axial)

            self.along_blade_quantities['force_x'].append(along_blade_force_x)
            self.along_blade_quantities['force_y'].append(along_blade_force_y)
            self.along_blade_quantities['force_z'].append(along_blade_force_z)

            # Create a form for quickly projecting turbine force onto a vector function space
            # equivalent to: tf_fn = project(tf, V, solver_type='cg')
            self.tf = tf
            self.tf_fn = Function(fs.V,name="tf_fn")

            # if isinstance(self.tf, list):
            #     self.blade_tf_fn = Function(fs.V)

            #     self.a = []
            #     self.L = []

            #     self.A = []
            #     self.b = []

            #     self.project_tf = []

            #     for blade_id in range(self.num_blades):
            #         blade_a = inner(TrialFunction(fs.V), TestFunction(fs.V))*dx
            #         blade_L = inner(self.tf[blade_id], TestFunction(fs.V))*dx

            #         self.a.append(blade_a)
            #         self.L.append(blade_L)

            #         self.A.append(assemble(self.a[blade_id]))
            #         self.b.append(assemble(self.L[blade_id]))

            #         self.project_tf.append(PETScKrylovSolver('cg', 'jacobi'))
            #         self.project_tf[blade_id].set_operator(self.A[blade_id])

            # else:
            #     self.a = inner(TrialFunction(fs.V), TestFunction(fs.V))*dx
            #     self.L = inner(self.tf, TestFunction(fs.V))*dx

            #     self.A = assemble(self.a)
            #     self.b = assemble(self.L)

            #     self.project_tf = PETScKrylovSolver('cg', 'jacobi')
            #     self.project_tf.set_operator(self.A)

            self.first_call_to_alm = False

        else:

            self.simTime_prev.assign(self.simTime)
            self.simTime.assign(kwargs['simTime'])

            self.platform_theta_prev.assign(self.platform_theta)
            platform_theta = self.calc_platform_motion()
            self.platform_theta.assign(platform_theta)

            # print(self.platform_theta_prev, self.platform_theta_prev.values())

            # TODO: we got to do something like this. By storing the aoa forms in a list, 
            # we should be able to just reassemble them after changing the time. We also 
            # might need to do something similar to velocity


            along_blade_lift = []
            along_blade_drag = []
            along_blade_aoa = []
            along_blade_axial = []

            along_blade_force_x = []
            along_blade_force_y = []
            along_blade_force_z = []

            for blade_id in range(self.num_blades):
                for actuator_id in range(self.num_actuator_nodes):
                    # Re Evaluate velocity
                    x_0_prev = self.x_0_prev[blade_id][actuator_id]

                    if self.use_gauss_vel_probe:
                        vel_fluid_temp = self.vel_fluid_eval(u, self.gauss_kernel[blade_id][actuator_id])
                    else:
                        vel_fluid_temp = mpi_eval(u, x_0_prev)

                    # vel_fluid_temp = [AdjFloat(8.0),AdjFloat(0.0),AdjFloat(0.0)]#mpi_eval(u, x_0_prev)
                    for i in range(self.dom.dim):
                        self.vel_fluid[blade_id][actuator_id][i].assign(vel_fluid_temp[i])

                    # Re Evaluate the aoa form
                    self.aoa_values[blade_id][actuator_id].assign(ufl_eval(self.aoa_forms[blade_id][actuator_id],print_statement=f"aoa: {(blade_id,actuator_id)}, time: {float(self.simTime):1.4f}, rank: {self.params.rank}"))

                    # Get blade props
                    rdim = self.rdim[actuator_id]
                    aoa = self.aoa_values[blade_id][actuator_id]

                    # Lookup the lift and drag coefficients
                    self.mcl[blade_id][actuator_id].assign(self.lookup_lift_coeff(rdim,aoa))
                    self.mcd[blade_id][actuator_id].assign(self.lookup_drag_coeff(rdim,aoa))

                    # print(f"lift_mag: {self.lift_mag[blade_id][actuator_id]}")
                    along_blade_lift.append(float(self.lift_mag[blade_id][actuator_id]))
                    along_blade_drag.append(float(self.drag_mag[blade_id][actuator_id]))
                    along_blade_aoa.append(float(self.aoa_values[blade_id][actuator_id]))
                    along_blade_axial.append(float(self.axial[blade_id][actuator_id]))

                    along_blade_force_x.append(self.actuator_force_components[blade_id][actuator_id][0])
                    along_blade_force_y.append(self.actuator_force_components[blade_id][actuator_id][1])
                    along_blade_force_z.append(self.actuator_force_components[blade_id][actuator_id][2])

            self.along_blade_quantities['lift'].append(along_blade_lift)
            self.along_blade_quantities['drag'].append(along_blade_drag)
            self.along_blade_quantities['aoa'].append(along_blade_aoa)
            self.along_blade_quantities['axial'].append(along_blade_axial)

            self.along_blade_quantities['force_x'].append(along_blade_force_x)
            self.along_blade_quantities['force_y'].append(along_blade_force_y)
            self.along_blade_quantities['force_z'].append(along_blade_force_z)

        # if hasattr(self, 'project_tf'):
        #     if isinstance(self.tf, list):
        #         self.fprint('Solving list w/ custom')
        #         self.tf_fn.vector()[:] = 0.0

        #         for blade_id in range(self.num_blades):
        #             self.b[blade_id] = assemble(self.L[blade_id], tensor=self.b[blade_id])
        #             self.project_tf[blade_id].solve(self.blade_tf_fn.vector(), self.b[blade_id])
        #             self.tf_fn.vector()[:] += self.blade_tf_fn.vector().get_local()
        #     else:
        #         self.fprint('Solving non-list w/ custom')
        #         self.b = assemble(self.L, tensor=self.b)
        #         self.project_tf.solve(self.tf_fn.vector(), self.b)

        # else:
        #     if isinstance(self.tf, list):
        #         self.fprint('Solving list w/ built in')

        #         self.tf_fn.vector()[:] = 0.0

        #         for blade_id in range(self.num_blades):
        #             self.blade_tf_fn = project(self.tf[blade_id], fs.V, solver_type='cg', preconditioner_type='jacobi')
        #             self.tf_fn.vector().axpy(1.0, self.blade_tf_fn.vector())

        #     else:
        #         self.fprint('Solving non-list w/ built in')
        #         self.tf_fn = project(self.tf, fs.V, solver_type='cg', preconditioner_type='jacobi')
        #         filename = os.path.join(self.params.folder, f'data/alm/true_force.npy')
        #         np.save(filename, self.tf_fn.vector()[:])

        # print(self.tf_fn.vector().max(), self.tf_fn.vector().min())

        # filename = os.path.join(self.params.folder, f'data/alm/true_force.npy')
        # truth_data = np.load(filename)

        # print(f'Max diff from truth: {np.amax(np.abs(truth_data - self.tf_fn.vector()[:]))}')

        # temp = 0
        # for blade_id in range(self.num_blades):
        #     temp += self.tf[blade_id]
        # self.tf_fn = project(temp, fs.V, solver_type='cg', preconditioner_type='jacobi')


        return sum(self.tf)
        # return self.tf_fn


    def power(self, u, inflow_angle):
        # Create a cylindrical expression aligned with the position of this turbine
        # self.cyld_expr = Expression(('sin(yaw)*(x[2]-zs)', '-cos(yaw)*(x[2]-zs)', '(x[1]-ys)*cos(yaw)-(x[0]-xs)*sin(yaw)'),
        #    degree=1,
        #    yaw=float(self.myaw),
        #    xs=float(self.mx),
        #    ys=float(self.my),
        #    zs=float(self.mz))

        #print(f'Theta = {float(self.platform_theta)}, mx = {float(self.mx)}, my = {float(self.my)}, mz = {float(self.mz)}')

        self.cyld_expr = Expression(
            ('-x[1]*sin(pitch)',
            '-1.0*(-sin(pitch)*x[0] + cos(pitch)*x[2] - zs)',
            'x[1]*cos(pitch)'),
            degree=1,
            pitch=-float(self.platform_theta),
            xs=float(self.mx),
            ys=float(self.my),
            zs=float(self.mz))

        # self.cyld_expr = Expression(('xs', '0.0', '0.0'), degree=1, xs=1.0)

        if isinstance(self.tf, list):
            self.power_dolfin_old = 0.0

            for blade_id in range(self.num_blades):
                self.power_dolfin_old += assemble(1e-6*dot(-self.tf[blade_id]*self.angular_velocity, self.cyld_expr)*dx)

        else:
            self.power_dolfin_old = assemble(1e-6*dot(-self.tf*self.angular_velocity, self.cyld_expr)*dx)

        # TODO: Assebmle these objects in groups to avoid running out of memory/hitting recursion limit
        # as a result of trying to assemble lots of actuator nodes at once.

        # rotor_torque = sum(tangential_force * rdim)
        # rotor_power = rotor_torque * (2*pi*RPM/60)
        rotor_torque = 0.0

        for blade_id in range(self.num_blades):
            for actuator_id in range(self.num_actuator_nodes):
                # torque = tangential_force * lever_arm
                rotor_torque += self.actuator_force_components[blade_id][actuator_id][2]*self.rdim[actuator_id]

        self.power_dolfin = assemble(-1e-6*rotor_torque*self.angular_velocity*dx)

        if self.params.rank == 0:
            print(f'old: {self.power_dolfin_old}')
            print(f'new: {self.power_dolfin}')

        # self.power_dolfin_old = assemble(dot(self.tf,as_vector((1.0,1.0,1.0)))*dx)
        # print('Power: ', self.power_dolfin_old)

        return self.power_dolfin

    def prepare_saved_functions(self, func_list):
        if len(func_list) == 0:
            func_list = [
                [self.tf_fn,"turbine_force"]
            ]
        else:
            func_list[1][0] += self.tf_fn

        return func_list

    def finalize_turbine(self):
        if self.params.rank == 0:
            for key in self.along_blade_quantities.keys():
                data = self.along_blade_quantities[key]
                try:
                    data = np.array(data)
                    filename = os.path.join(self.params.folder, f'data/alm/{key}_{self.index}.npy')
                    np.save(filename, data)
                except:
                    self.fprint(f'Could not save along blade data: {key}')
