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


    def init_constant_alm_terms(self, fs):

        #================================================================
        # Get information about the mesh coordinates and distances
        #================================================================

        self.ndim = self.dom.dim

        # Get the coordinates of the vector function space
        self.coords_base = fs.V.tabulate_dof_coordinates()
        self.coords_base = np.copy(self.coords_base[0::self.ndim, :])

        # Resape a linear copy of the coordinates for every mesh point
        self.coords = np.copy(self.coords_base)
        self.coordsLinear = np.copy(self.coords.reshape(-1, 1))

        bbox = self.dom.mesh.bounding_box_tree()
        turbine_loc_point = Point(self.x, self.y, self.z)
        node_id, dist = bbox.compute_closest_entity(turbine_loc_point)
        self.min_dist = dist

        #================================================================
        # Set turbine properties and calculate derived values
        #================================================================

        # Set the number of blades in the turbine
        self.num_blades = 3

        # Set the spacing pf each blade
        self.theta_0_vec = np.linspace(0.0, 2.0*np.pi, self.num_blades, endpoint = False)

        # Calculate the blade velocity
        self.angular_velocity = 2.0*np.pi*self.rpm/60.0
        self.tip_speed = self.angular_velocity*self.radius

        # Recommendation from Churchfield et al.
        # self.gaussian_width = 2.0*0.035*2.0*self.radius
        self.gaussian_width = float(self.gauss_factor)*self.dom.global_hmin

        if self.blade_segments == "computed":
            self.num_blade_segments = int(2.0*self.radius/self.gaussian_width)
        else:
            self.num_blade_segments = self.blade_segments

        # Calculate the radial position of each actuator node
        self.rdim = np.linspace(0.0, self.radius, self.num_blade_segments)

        # Calculate width of an individual blade segment
        # w = rdim[1] - rdim[0]
        self.w = (self.rdim[1] - self.rdim[0])*np.ones(self.num_blade_segments)
        self.w[0] = self.w[0]/2.0
        self.w[-1] = self.w[-1]/2.0

        #================================================================
        # Set constant structures
        #================================================================

        # Calculate an array describing the x, y, z position of each actuator node
        # Note: The basic blade is oriented along the +y-axis
        self.blade_pos_base = np.vstack((np.zeros(self.num_blade_segments),
                                         self.rdim,
                                         np.zeros(self.num_blade_segments)))

        # Specify the velocity vector at each actuator node
        # Note: A blade with span oriented along the +y-axis moves in the +z direction
        # Note: blade_vel_base is negative since we seek the velocity of the fluid relative to a stationary blade
        # and blade_vel_base is defined based on the movement of the blade
        self.blade_vel_base = np.vstack((np.zeros(self.num_blade_segments),
                                         np.zeros(self.num_blade_segments),
                                         np.linspace(0.0, -self.tip_speed, self.num_blade_segments)))

        # Create unit vectors aligned with blade geometry
        # blade_unit_vec_base[:, 0] = points along rotor shaft
        # blade_unit_vec_base[:, 1] = points along blade span axis
        # blade_unit_vec_base[:, 2] = points tangential to blade span axis (generates a torque about rotor shaft)
        self.blade_unit_vec_base = np.array([[1.0, 0.0, 0.0],
                                             [0.0, 1.0, 0.0],
                                             [0.0, 0.0, 1.0]])

        #================================================================
        # Finally, initialize important dolfin terms
        #================================================================

        # Create a Constant "wrapper" to enable dolfin to track mpi_u_fluid
        # self.mpi_u_fluid_constant = Constant(np.zeros(self.ndim*self.num_blades*self.num_blade_segments), name="mpi_u_fluid")


        # self.mchord.append(turb_i_chord)
        # self.mtwist.append(turb_i_twist)
        # self.mcl.append(turb_i_lift)
        # self.mcd.append(turb_i_drag)
        # self.chord = np.array(self.mchord,dtype=float)
        # self.cl = np.array(self.mcl,dtype=float)
        # self.cd = np.array(self.mcd,dtype=float)
        # self.farm.baseline_chord = np.array(self.chord[0])/self.chord_factor

        # self.cyld_expr_list = [None]*self.farm.numturbs

        # FIXME: need to get these coordinates the correct way
        # Make this a list of constants

        # self.mchord = Constant(self.chord, name="chord_{:d}".format(self.index))
        # self.mtwist = Constant(self.twist, name="twist_{:d}".format(self.index))
        # self.mlift = Constant(self.lift, name="lift_{:d}".format(self.index))
        # self.mdrag = Constant(self.drag, name="drag_{:d}".format(self.index))


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
            interp_points = np.linspace(0.0, 1.0, self.num_blade_segments)

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
            chord = self.radius/20.0*np.ones(self.num_blade_segments)
            twist = np.zeros(self.num_blade_segments)
            cl = np.ones(self.num_blade_segments)
            cd = 0.1*np.ones(self.num_blade_segments)

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

                s_0 = station_radii[station_id]*np.ones(np.shape(data)[0])

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
                    station = np.hstack((station, s_i))
                    angles = np.hstack((angles, angles_i))
                    c_lift = np.hstack((c_lift, c_lift_i))
                    c_drag = np.hstack((c_drag, c_drag_i))

            nodes = np.vstack((station, angles)).T

            # Create interpolation functions for lift and drag based on angle of attack and location along blade
            self.lookup_cl = interp.LinearNDInterpolator(nodes, c_lift)
            self.lookup_cd = interp.LinearNDInterpolator(nodes, c_drag)

        else:
            self.lookup_cl = None
            self.lookup_cd = None


    def init_unsteady_alm_terms(self):

        '''
        
        This function does the resetting of summation-like terms for
        calculating the rotor and handles the time stepping for where the 
        actuator forces are applied and where the fluid velocity is probed.

        simTime_prev: the time value associated with the last fluid solve
        simTime: the time value associated with the current fluid solve

        theta_prev: the rotation of the rotor at time 0.5*(simTime_prev + simTime),
        i.e., halfway between simTime_prev and simTime
        theta: the rotation of the rotor at time (simTime + 0.5*dt),
        i.e., halfway between the current simTime and next simTime

        simTime_prev                   simTime             simTime+dt
              |        theta_prev         |        theta        |
              |             |             |          |          |
              |             |             |          |          |
              ---------------------------------------------------
              |<--------- dt_prev ------->|<-------- dt ------->|
              |                           |                     |
             k-1                          k                    k+1

        '''

        self.rotor_torque = np.zeros(1)
        self.rotor_torque_count = np.zeros(1, dtype=int)

        # The fluid velocity should be probed at the current position/time MINUS 0.5*dt
        if self.simTime_prev is None:
            self.theta_prev = self.simTime*self.angular_velocity
        else:
            self.theta_prev = 0.5*(self.simTime_prev + self.simTime)*self.angular_velocity

        self.theta = (self.simTime + 0.5*self.dt)*self.angular_velocity




    def update_alm_node_positions(self):

        def wave_motion(tt):

            # Peak-to-peak duration of one heave cycle (seconds)
            w = 10.0

            # Maximum amplitude of a heave cycle (degrees)
            # a = 20.0
            a = 0.0

            wave_theta = np.sin(2.0*np.pi*tt/w)
            wave_theta *= np.radians(a)

            return wave_theta

        self.blade_pos = []
        self.blade_pos_prev = []
        self.blade_unit_vec = []
        self.blade_vel = []

        init_points = [self.blade_pos_base, self.blade_unit_vec_base, self.blade_vel_base]

        for blade_id, theta_0 in enumerate(self.theta_0_vec):
            # Rotate the blade into the correct position around the x-axis (due to rotor spin)
            pos, unit_vec, vel = self.rotate_points([self.blade_pos_base, self.blade_unit_vec_base, self.blade_vel_base], 
                                                     self.theta + theta_0, [1, 0, 0])

            pos_prev = self.rotate_points(self.blade_pos_base, self.theta_prev + theta_0, [1, 0, 0])


            # Rotate the blade into the correct position around the z-axis (due to yaw)
            pos, pos_prev, unit_vec, vel = self.rotate_points([pos, pos_prev, unit_vec, vel], float(self.myaw), [0, 0, 1])


            # Translate the positions into the correct location
            pos, pos_prev = self.translate_points([pos, pos_prev], [self.x, self.y, self.z])


            # Rotate the points due to the waves
            theta_tt = self.theta/self.angular_velocity
            theta_prev_tt = self.theta_prev/self.angular_velocity

            pos, unit_vec, vel = self.rotate_points([pos, unit_vec, vel], wave_motion(theta_tt), [0, 1, 0])

            pos_prev = self.rotate_points(pos_prev, wave_motion(theta_prev_tt), [0, 1, 0])

            # if self.turbine_motion_freq is not None:
            #     motion_theta = self.turbine_motion_amp*np.sin(self.turbine_motion_freq*self.simTime_ahead*np.pi*2.0)
            #     motion_theta = self.turbine_motion_amp*np.sin(self.turbine_motion_freq*self.simTime_behind*np.pi*2.0)
            #     Ry = self.rotation_mat_about_axis(motion_theta, [0, 1, 0])
            #     blade_pos = np.dot(Ry, blade_pos)
            #     self.blade_pos_prev.append(np.copy(blade_pos))
            # if self.turbine_motion_freq is not None:
            #     motion_theta = self.turbine_motion_amp*np.sin(self.turbine_motion_freq*self.simTime_ahead*np.pi*2.0)
            #     Ry = self.rotation_mat_about_axis(motion_theta, [0, 1, 0])
            #     self.blade_vel[blade_id] = np.dot(Ry, self.blade_vel[blade_id])
            #     self.blade_unit_vec[blade_id] = np.dot(Ry, self.blade_unit_vec[blade_id])
            #     self.blade_pos[blade_id] = np.dot(Ry, self.blade_pos[blade_id])

            #     turbine_motion_vel = (self.blade_pos[blade_id] - self.blade_pos_prev[blade_id])/(self.simTime_ahead - self.simTime_behind)

            #     if blade_id == 0:
            #         print('Hub Location:', blade_pos[:, 0])
            #         print('Theta:', motion_theta)

            self.blade_pos.append(pos)
            self.blade_unit_vec.append(unit_vec)
            self.blade_vel.append(vel)
            self.blade_pos_prev.append(pos_prev)

    def get_u_fluid_at_alm_nodes(self, u_k):

        # Create an empty list to hold all the components of velocity
        mpi_u_fluid = []

        for blade_id in range(self.num_blades):
            # Need to probe the velocity point at each actuator node,
            # where actuator nodes are individual columns of blade_pos

            blade_mpi_u_fluid = []

            for k in range(self.num_blade_segments):
                pos_i = np.array([0.0, 0.0, 0.0])

                # If using the local velocity, measure at the blade
                if self.use_local_velocity:
                    pos_i[0] = self.blade_pos_prev[blade_id][0, k]
                else:
                    pos_i[0] = self.dom.x_range[0]

                pos_i[1] = self.blade_pos_prev[blade_id][1, k]
                pos_i[2] = self.blade_pos_prev[blade_id][2, k]

                # Try to access the fluid velocity at this actuator point
                # If this rank doesn't own that point, an error will occur,
                # in which case zeros should be reported
                blade_mpi_u_fluid.append(mpi_eval(u_k, pos_i))

            mpi_u_fluid.append(blade_mpi_u_fluid)

        # Store the finished list of Constants
        self.mpi_u_fluid = mpi_u_fluid


    def lookup_lift_and_drag(self, u_rel, twist, unit_vec):

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

        # Initialize the real cl and cd profiles
        real_cl = np.zeros(self.num_blade_segments)
        real_cd = np.zeros(self.num_blade_segments)

        # fp = open(problem.aoa_file, 'a')

        tip_loss_fac = np.zeros(self.num_blade_segments)

        aoa_list = []

        for k in range(self.num_blade_segments):
            # Get the relative wind velocity at this node
            wind_vec = u_rel[:, k]

            # Remove the component in the radial direction (along the blade span)
            wind_vec -= np.dot(wind_vec, unit_vec[:, 1])*unit_vec[:, 1]

            # aoa = get_angle_between_vectors(arg1, arg2, arg3)
            # arg1 = in-plane vector pointing opposite rotation (blade sweep direction)
            # arg2 = relative wind vector at node k, including blade rotation effects (wind direction)
            # arg3 = unit vector normal to plane of rotation, in this case, radially along span
            aoa = get_angle_between_vectors(-unit_vec[:, 2], wind_vec, -unit_vec[:, 1])

            # Compute tip-loss factor
            if self.tip_loss:
                if self.rdim[k] < 1e-12:
                    tip_loss_fac[k] = 1.0

                else:
                    loss_exponent = 3.0/2.0*(self.radius-self.rdim[k])/(self.rdim[k]*np.sin(aoa))
                    acos_arg = np.exp(-loss_exponent)
                    acos_arg = np.clip(acos_arg, -1.0, 1.0)
                    tip_loss_fac[k] = 2.0/np.pi*np.arccos(acos_arg)

            else:
                tip_loss_fac[k] = 1.0

            # Remove the portion of the angle due to twist
            aoa -= twist[k]

            if self.hub_rad > 0.0:
                if self.rdim[k] - self.hub_rad < 0:
                    lookup_rdim = 0.0
                else:
                    lookup_rdim = self.radius*(self.rdim[k] - self.hub_rad)/(self.radius - self.hub_rad)

            else:
                lookup_rdim = self.rdim[k]

            # Store the cl and cd by interpolating this (aoa, span) pair from the tables
            real_cl[k] = self.lookup_cl(lookup_rdim, aoa)
            real_cd[k] = self.lookup_cd(lookup_rdim, aoa)

            aoa_list.append(aoa)

            # Write the aoa to a file for future reference
            # fa.write('%.5f, ' % (aoa/np.pi*180.0))

        aoa_list = np.array(aoa_list)

        # fp.close()
        return real_cl, real_cd, tip_loss_fac, aoa_list


    def build_actuator_node(self, u_k, x_0, orientation):

        gauss_function = Function(self.fs.V)

        x = SpatialCoordinate(self.dom.mesh)

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
            self._setup_chord_force()
            chord = self.mchord[i]
            force = self.radial_chord_force(r,chord)
        F = -0.5*A*C_tprime*force

        ### Calculate normalization constant ###
        volNormalization = T_norm*D_norm*W*R**(self.dom.dim-1)


        return gauss_function
 


    def finalize_rotor_torque(self):

        data_in_torque = np.zeros(self.params.num_procs)
        self.params.comm.Gather(self.rotor_torque, data_in_torque, root=0)

        data_in_torque_count = np.zeros(self.params.num_procs, dtype=int)
        self.params.comm.Gather(self.rotor_torque_count, data_in_torque_count, root=0)

        if self.params.rank == 0:
            rotor_torque_sum = np.sum(data_in_torque)
            rotor_torque_count_sum = np.sum(data_in_torque_count)

            # This removes the possibility of a power being doubled or tripled
            # if multiple ranks include this turbine and therefore calculate a torque
            self.rotor_torque = rotor_torque_sum/rotor_torque_count_sum

        self.params.comm.Bcast(self.rotor_torque, root=0)


    def turbine_force(self, u, inflow_angle, fs, **kwargs):
        # if dfd is None, alm_output is a dolfin function (tf) [1 x numPts*ndim]
        # otherwise, it returns a numpy array of derivatives [numPts*ndim x numControls]

        try:
            self.simTime = kwargs['simTime']
            self.simTime_prev = kwargs['simTime_prev']
            self.dt = kwargs['dt']
        except:
            raise ValueError('"simTime" and "dt" must be specified for the calculation of ALM force.')

        # If this is the first call to the function, set some things up before proceeding
        if self.first_call_to_alm:
            self.init_constant_alm_terms(fs)
            self.init_blade_properties()
            self.init_lift_drag_lookup_tables()
            self.create_controls(initial_call_from_setup=False)

            self.fprint(f'Gaussian Width: {self.gaussian_width}')
            self.fprint(f'Minimum Mesh Spacing: {self.dom.global_hmin}')
            self.fprint(f'Number of Actuators per Blade: {self.num_blade_segments}')

        # Initialize summation, counting, etc., variables for alm solve
        self.init_unsteady_alm_terms()

        # self.global_to_local_coord_transform()
        self.update_alm_node_positions()

        # Call the function to build the complete mpi_u_fluid array
        self.get_u_fluid_at_alm_nodes(u)

        # Call the ALM function for this turbine
        self.tf = self.build_actuator_lines(u, inflow_angle, fs)

        # Do some sharing of information when everything is finished
        self.finalize_rotor_torque()

        if self.first_call_to_alm:
            self.first_call_to_alm = False

        return self.tf


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
























from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
# from windse.helper_functions import mpi_eval
# import scipy.interpolate as interp

def fake_scipy_lift(rdim, aoa):
    
    return 1.0


def fake_scipy_drag(rdim, aoa):
    
    return 0.1


def mpi_eval(u_k, x_0):
    a = Constant(u_k(x_0))
    
    return a


def rot_x(theta):
    Rx = as_tensor([[1, 0, 0],
                   [0, cos(theta), -sin(theta)],
                   [0, sin(theta), cos(theta)]])
    
    return Rx


def calculate_relative_fluid_velocity(u_k, x_0, n_0, rdim):

    # Evaluate the fluid velocity u_k at the point x_0 
    vel_fluid = mpi_eval(u_k, x_0)
    
    # Remove the component of the fluid velocity oriented along the blade's axis
    vel_fluid -= dot(vel_fluid, n_0[:, 1])*n_0[:, 1]

    # TODO : replace with angular velocity calc
    angular_velocity = 9.0
    vel_blade = rdim*angular_velocity*n_0[:, 2]

    # The relative velocity is the sum of the fluid and blade velocity
    vel_rel = vel_fluid + vel_blade
    
    # Calculate the magnitude of the relative velocity
    vel_rel_mag = sqrt(vel_rel[0]**2.0 + vel_rel[1]**2.0 + vel_rel[2]**2.0)
    vel_rel_unit = vel_rel/vel_rel_mag

    return vel_rel, vel_rel_mag, vel_rel_unit


def calcuate_aoa(vel_rel, n_0):
        
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


def lookup_lift_coeff(rdim, aoa):
    
#     TODO: Jeff, make this take in rdim and aoa, analogous to looking up terrain 
#     val = self.dom.Ground(float(x),float(y),dx=dx,dy=dy)
#     return Constant(val)
    cl = fake_scipy_lift(rdim, aoa)

    return Constant(cl)


def lookup_drag_coeff(rdim, aoa):
    
#     TODO: Jeff, make this take in rdim and aoa, analogous to looking up terrain 
#     val = self.dom.Ground(float(x),float(y),dx=dx,dy=dy)
#     return Constant(val)
    cd = fake_scipy_drag(rdim, aoa)
    
    return Constant(cd)


def build_lift_drag_vec(vel_rel_unit, n_0):
    # The drag unit vector is oriented opposite of the relative fluid velocity
    drag_unit = -vel_rel_unit
    
    # The lift unit vector is perpendicular to the drag unit vector and the blade axis
    lift_unit = cross(drag_unit, n_0[:, 1])
    
    return lift_unit, drag_unit


def build_gauss_kernel_at_point(u_k, x_0, eps):

    # Get the coordinates of the mesh points as UFL term
    x_mesh = SpatialCoordinate(u_k.function_space().mesh())
    
    # Calculate distance from actuator node, x_0, in 3 directions
    delta_x = x_mesh[0] - x_0[0]
    delta_y = x_mesh[1] - x_0[1]
    delta_z = x_mesh[2] - x_0[2]
    
    # Calculate the distance squared, i.e., (|x_mesh-x|)^2
    r2 = delta_x**2.0 + delta_y**2.0 + delta_z**2.0
    
    # Compute the Gaussian kernel
    gauss_kernel = exp(-r2/(eps**2.0))/(eps**3.0*np.pi**1.5)

    return gauss_kernel


def calculate_tip_loss(rdim, aoa):

    tip_loss_fac = 1.0

    return tip_loss_fac


def build_actuator_node(u_k, x_0, n_0, rdim):
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
    vel_rel, vel_rel_mag, vel_rel_unit = calculate_relative_fluid_velocity(u_k, x_0, n_0, rdim)

    # Calculate the angle of attack for this actuator node node
    aoa = calcuate_aoa(vel_rel, n_0)
    
    # Lookup the lift and drag coefficients
    cl = lookup_lift_coeff(rdim, aoa)
    cd = lookup_drag_coeff(rdim, aoa)

    # Calculate the tip loss factor
    tip_loss_fac = calculate_tip_loss(rdim, aoa)
            
    # Calculate the magnitude (a scalar quantity) of the lift and drag forces
    lift_mag = tip_loss_fac*(0.5*cl*rho*chord*w*vel_rel_mag**2.0)
    drag_mag = tip_loss_fac*(0.5*cd*rho*chord*w*vel_rel_mag**2.0)

    # Build a spherical Gaussian with width eps about point x_0
    gauss_kernel = build_gauss_kernel_at_point(u_k, x_0, eps)
    
    # Smear the lift and drag forces onto the grid using the Gaussian kernel
    # The result is a scalar function at each point on the grid
    projected_lift_mag = lift_mag*gauss_kernel
    projected_drag_mag = drag_mag*gauss_kernel
    
    # Compute the orientation of the lift and drag unit vectors
    lift_unit, drag_unit = build_lift_drag_vec(vel_rel_unit, n_0)
    
    # The smeared scalar force magnitude is oriented along the lift and drag directions
    lift_force = projected_lift_mag*lift_unit
    drag_force = projected_drag_mag*drag_unit

    # The final force of this actuator is the sum of lift and drag effects
    actuator_force = lift_force + drag_force
        
    return actuator_force


mesh = UnitCubeMesh(10, 10, 10)

V = VectorFunctionSpace(mesh, 'P', 1)

u_k = Function(V)

e = Expression(('x[0]', 'x[1]', 'x[2]'), degree=1)

u_k.interpolate(e)

# This represents the actual evaluation point
# (after shifting by hub height, location, yaw, platform mvmt, etc.)
x_0 = as_vector([0.5, 0.6, 0.7])

# This represents the actual orientation of the blade as a unit vector
# (after rotating, yawing, floating, etc.)
n_0 = as_tensor([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]])

# # Orient this by dotting with lift direction
# rx = rot_x(np.radians(90))

chord = 1.0
rho = 1.0
w = 1.0
eps = 0.2

rdim = 0.5

# for all actuator nodes:
print('starting actuator node')
g = build_actuator_node(u_k, x_0, n_0, rdim)
print('finished')


