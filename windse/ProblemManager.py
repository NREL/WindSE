"""
The ProblemManager contains all of the 
different classes of problems that windse can solve
"""

import __main__
import os

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np
    import time
    import scipy.interpolate as interp
    import glob

    ### Import the cumulative parameters ###
    from windse import windse_parameters
    # from memory_profiler import memory_usage

    ### Check if we need dolfin_adjoint ###
    if windse_parameters.dolfin_adjoint:
        from dolfin_adjoint import *

class GenericProblem(object):
    """
    A GenericProblem contains on the basic functions required by all problem objects.
    
    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_data):
        ### save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.dom  = domain
        self.farm = windfarm
        self.fs   = function_space 
        self.bd  = boundary_data
        self.tf_first_save = True
        self.fprint = self.params.fprint
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode

        ### Update attributes based on params file ###
        for key, value in self.params["problem"].items():
            setattr(self,key,value)

        self.record_time = self.params["optimization"].get("record_time",0.0)
        if isinstance(self.record_time,str):
            self.record_time = 0.0

    def DebugOutput(self):
        if self.debug_mode:

            # integral of nu_t
            int_nut = assemble(self.nu_T*dx)/self.dom.volume
            self.tag_output("int_nu_T", int_nut)

            # integral of tf
            if self.dom.dim == 3:
                e1 = Constant((1,0,0)); e2 = Constant((0,1,0)); e3 = Constant((0,0,1));
            else:
                e1 = Constant((1,0)); e2 = Constant((0,1));

            int_tf_x = assemble(inner(self.tf,e1)*dx)/self.dom.volume
            self.tag_output("int_tf_x", int_tf_x)
            int_tf_y = assemble(inner(self.tf,e2)*dx)/self.dom.volume
            self.tag_output("int_tf_y", int_tf_y)
            if self.dom.dim == 3:
                int_tf_z = assemble(inner(self.tf,e3)*dx)/self.dom.volume
                self.tag_output("int_tf_z", int_tf_z)


            if self.farm.turbine_method == 'alm':
                self.tag_output("min_chord", np.min(self.chord))
                self.tag_output("max_chord", np.max(self.chord))
                self.tag_output("avg_chord", np.mean(self.chord))
                self.tag_output("min_cl", np.min(self.cl))
                self.tag_output("max_cl", np.max(self.cl))
                self.tag_output("avg_cl", np.mean(self.cl))
                self.tag_output("min_cd", np.min(self.cd))
                self.tag_output("max_cd", np.max(self.cd))
                self.tag_output("avg_cd", np.mean(self.cd))
                self.tag_output("num_blade_segments", self.num_blade_segments)

    def ComputeTurbineForce(self,u,inflow_angle,simTime=0.0):

        ### Compute the relative yaw angle ###
        if inflow_angle is not None:
            inflow_angle = inflow_angle-self.dom.inflow_angle
        else:
            inflow_angle = self.dom.inflow_angle

        self.fprint('Computing turbine forces using %s' % (self.farm.turbine_method.upper()))

        ### Create the turbine force function ###
        if self.farm.turbine_method == "disabled" or self.farm.numturbs == 0:
            tf = Function(self.fs.V)
            self.farm.actuator_disks = []

        elif self.farm.turbine_method == "dolfin":
            self.tf1, self.tf2, self.tf3 = self.farm.DolfinTurbineForce(self.fs,self.dom.mesh,inflow_angle=inflow_angle)
            self.num_blade_segments = self.farm.blade_segments
            tf = self.tf1*u[0]**2+self.tf2*u[1]**2+self.tf3*u[0]*u[1]

        elif self.farm.turbine_method == "numpy":
            self.tf1, self.tf2, self.tf3 = self.farm.NumpyTurbineForce(self.fs,self.dom.mesh,inflow_angle=inflow_angle)
            tf = -(self.tf1*u[0]**2+self.tf2*u[1]**2+self.tf3*u[0]*u[1]) #negative or otherwise we get jets

        elif self.farm.turbine_method == 'alm':
            self.rpm = self.params["wind_farm"]["rpm"]

            hmin = self.dom.mesh.hmin()/np.sqrt(3)
            # self.gaussian_width = 2.0*hmin/3.0
            # Recommendation from Churchfield et al.
            self.gaussian_width = 2.0*0.035*2.0*self.farm.radius[0]
            self.gaussian_width = float(self.params["wind_farm"]["gauss_factor"])*hmin

            self.chord_factor = float(self.params["wind_farm"]["chord_factor"])

            self.fprint('Minimum Space Between Mesh: %f' % (hmin))
            self.fprint('Gaussian Width: %f' % (self.gaussian_width))

            # self.num_blade_segments = 10
            # self.num_blade_segments = int(10.0*self.farm.radius[0]/hmin)
            if self.farm.blade_segments == "computed":
                self.num_blade_segments = int(2.0*self.farm.radius[0]/self.gaussian_width)
                self.farm.blade_segments = self.num_blade_segments
            else:
                self.num_blade_segments = self.farm.blade_segments
            self.fprint('Num blade segments: %d' % (self.num_blade_segments))

            self.mchord = []
            self.mtwist = []
            self.mcl = []
            self.mcd = []

            self.num_times_called = 0
            self.first_call_to_alm = True
            self.blade_pos_previous = [[], [], []]
            self.simTime_list = []
            self.dt_list = []
            self.rotor_torque_dolfin_time = []
            self.simTime_id = 0

            ### create output files for alm data ###
            force_folder = self.params.folder+"data/alm/rotor_force/"
            aoa_folder = self.params.folder+"data/alm/angle_of_attack/"
            if not os.path.exists(aoa_folder): os.makedirs(aoa_folder)
            if not os.path.exists(force_folder): os.makedirs(force_folder)
            self.aoa_files = []
            self.force_files = []
            for i in range(self.farm.numturbs):
                self.aoa_files.append(aoa_folder+"aoa_turb_"+repr(i)+".csv")
                temp = ["x","y","z"]
                self.force_files.append([])
                for j in range(3):
                    self.force_files[i].append(force_folder+temp[j]+"_force_turb_"+repr(i)+".csv")

                ### Create header for angle of attack 
                fp = open(self.aoa_files[i], 'w')
                fp.write('time, ')
                for j in range(3):
                    for k in range(self.num_blade_segments):
                        fp.write('r%d_n%03d, ' % (j, k))
                fp.write('\n')
                fp.close()

                ### Create header for force files ###
                for dir_file in self.force_files[i]:
                    fp = open(dir_file, 'w')
                    fp.write('time, ')
                    for j in range(3):
                        for k in range(self.num_blade_segments):
                            fp.write('r%d_n%03d, ' % (j, k))
                    fp.write('\n')
                    fp.close()

            # Initialize the lift and drag files
            for fn in ['lift', 'drag']:
                fp = open('./output/%s/nodal_%s.csv' % (self.params.name, fn), 'w')
                fp.write('sim_time, theta, ')

                for j in range(self.num_blade_segments):
                    if j < self.num_blade_segments-1:
                        fp.write('%s_%02d, ' % (fn, j))
                    else:
                        fp.write('%s_%02d\n' % (fn, j))
                fp.close()

            turb_data = self.params["wind_farm"]["read_turb_data"]

            if turb_data:

                def build_lift_and_drag_tables(airfoil_data_path):

                    # Determine the number of files in airfoil_data_path
                    num_files = len(glob.glob('%s/*.txt' % (airfoil_data_path)))

                    for file_id in range(num_files):
                        # print('Reading Airfoil Data #%d' % (file_id))
                        data = np.genfromtxt('%s/af_station_%d.txt' % (airfoil_data_path, file_id), skip_header=1, delimiter=' ')

                        if file_id == 0:
                            # If this is the first file, store the angle data
                            interp_angles = data[:, 0]
                            num_angles = np.size(interp_angles)
                            
                            # If this is the first file, allocate space for the tables        
                            lift_table = np.zeros((num_angles, num_files))
                            drag_table = np.zeros((num_angles, num_files))
                            
                        # Store all the lift and drag data in the file_id column
                        lift_table[:, file_id] = data[:, 1]
                        drag_table[:, file_id] = data[:, 2]

                    return lift_table, drag_table, interp_angles


                self.fprint('Setting chord, lift, and drag from file \'%s\'' % (turb_data))

                actual_turbine_data = np.genfromtxt(turb_data, delimiter = ',', skip_header = 1)

                actual_x = actual_turbine_data[:, 0]

                actual_chord = self.chord_factor*actual_turbine_data[:, 1]

                # Baseline twist is expressed in degrees, convert to radians
                actual_twist = actual_turbine_data[:, 2]/180.0*np.pi

                actual_cl = actual_turbine_data[:, 3]
                actual_cd = actual_turbine_data[:, 4]


                lift_table, drag_table, interp_angles = build_lift_and_drag_tables('airfoil_polars')
                self.lift_table = lift_table
                self.drag_table = drag_table
                self.interp_angles = interp_angles


                modify_chord = False

                if modify_chord:
                    shift_amt = 0.25
                    low_end = 1.0 - shift_amt
                    high_end = 1.0 + shift_amt
                    chord_shift_amt = np.linspace(low_end, high_end, np.size(actual_chord))
                    actual_chord = chord_shift_amt*actual_chord

                # print('chord measured: ', actual_chord)
                # print('lift measured: ', actual_cl)
                # print('drag measured: ', actual_cd)

                # Create interpolators for chord, lift, and drag
                chord_interp = interp.interp1d(actual_x, actual_chord)
                twist_interp = interp.interp1d(actual_x, actual_twist)
                cl_interp = interp.interp1d(actual_x, actual_cl)
                cd_interp = interp.interp1d(actual_x, actual_cd)

                if self.params['problem']['script_iterator'] == 1:
                    # actual_chord_override = np.array([4.749999999999964473e+00,
                    #     4.749999999999989342e+00,
                    #     4.750000000000000000e+00,
                    #     4.749999999999996447e+00,
                    #     4.749999999999998224e+00,
                    #     4.749999999999999112e+00,
                    #     4.450663333333332083e+00,
                    #     3.867691040375480060e+00,
                    #     8.488750327288195896e-01,
                    #     1.000000000000000056e-01])

                    # Rollercoaster chord
                    # actual_chord_override = np.array([2.625003373758692859e+00,
                    #     4.155334586096127047e+00,
                    #     5.350334971114000915e+00,
                    #     5.349201385073955706e+00,
                    #     4.113615408074356772e+00,
                    #     2.310353854848651967e+00,
                    #     1.275891949427846672e+00,
                    #     1.360607260649611794e+00,
                    #     2.500724540479074776e+00,
                    #     2.000000000000000111e-01])

                    # two-turbine power optimization result
                    # actual_chord_override = np.array([2.599932521058577972e+00,
                    #     3.655152166735353259e+00,
                    #     4.303587816194747617e+00,
                    #     3.987666943125669050e+00,
                    #     3.191648166607643500e+00,
                    #     2.553086977597530449e+00,
                    #     2.242893315052848724e+00,
                    #     2.172796739396936960e+00,
                    #     2.000000000000000111e-01])

                    # two-turbine power optimization result (newer result, obtained using cfl = 2.0 instead of 2.5)
                    actual_chord_override = np.array([2.599908825,
                        3.654252701,
                        4.300236963,
                        3.978497397,
                        3.178777094,
                        2.545014297,
                        2.247624970,
                        2.188315504,
                        0.200000000])

                    actual_x_override = np.linspace(0.0, 1.0, np.size(actual_chord_override))

                    chord_interp_override = interp.interp1d(actual_x_override, actual_chord_override)

                if self.params['problem']['script_iterator'] == 2:
                    # actual_chord_override = np.array([5.200000000000000178e+00,
                    #     7.027575555710222410e+00,
                    #     8.476945227498621449e+00,
                    #     8.280660000000001020e+00,
                    #     6.990223354798463795e+00,
                    #     5.499632522453062222e+00,
                    #     4.450663333333333860e+00,
                    #     3.867691040375480060e+00,
                    #     8.488750327288112629e-01,
                    #     1.000000000000000056e-01])

                    # two-turbine power optimization result
                    # actual_chord_override = np.array([2.599932521058577972e+00,
                    #     3.655152166735353259e+00,
                    #     4.303587816194747617e+00,
                    #     3.987666943125669050e+00,
                    #     3.191648166607643500e+00,
                    #     2.553086977597530449e+00,
                    #     2.242893315052848724e+00,
                    #     2.172796739396936960e+00,
                    #     2.000000000000000111e-01])

                    # actual_chord_override[1:-1] = 1.05*actual_chord_override[1:-1]


                    # two-turbine power optimization result (with 50% more change along optimized direction)
                    # actual_chord_override = np.array([2.59989878,
                    #     3.656438255,
                    #     4.312742105,
                    #     4.00897182,
                    #     3.235834005,
                    #     2.632667705,
                    #     2.37238998,
                    #     2.381786645,
                    #     0.2])

                    # two-turbine power optimization result (with 50% more change along optimized direction)
                    # newer result, obtained using cfl = 2.0 instead of 2.5
                    actual_chord_override = np.array([2.599863238,
                        3.655089051,
                        4.307715820,
                        3.995217504,
                        3.216527391,
                        2.620558679,
                        2.379487455,
                        2.405064793,
                        0.200000000])

                    actual_x_override = np.linspace(0.0, 1.0, np.size(actual_chord_override))

                    chord_interp_override = interp.interp1d(actual_x_override, actual_chord_override)

                # Construct the points at which to generate interpolated values
                interp_points = np.linspace(0.0, 1.0, self.num_blade_segments)

                # Generate the interpolated values
                chord = chord_interp(interp_points)
                twist = twist_interp(interp_points)
                cl = cl_interp(interp_points)
                cd = cd_interp(interp_points)

                if self.params['problem']['script_iterator'] > 0:
                    chord_override = chord_interp_override(interp_points)

            else:
                # If not reading from a file, prescribe dummy values
                chord = self.farm.radius[0]/20.0*np.ones(self.num_blade_segments)
                twist = np.zeros(self.num_blade_segments)
                cl = np.ones(self.num_blade_segments)
                cd = 0.1*np.ones(self.num_blade_segments)

            # Save the list of controls to self
            for turb_i in range(self.farm.numturbs):
                turb_i_chord = []
                turb_i_twist = []
                turb_i_lift = []
                turb_i_drag = []

                if turb_i == 0 and self.params['problem']['script_iterator'] > 0:
                    chord_to_use = chord_override
                else:
                    # chord_to_use = [2.5993912952866394, 3.5492808721592186, 4.67042880238879, 4.94984740731984, 4.791894547336281, 4.520032008634496, 4.450663333333334, 3.86769104037548, 3.395500130915245, 0.2]
                    # chord_to_use = [2.59978608485722,    3.53358366863681,    4.44802950633928,    4.50303965573438,    4.23785738559839,    4.06900117097451,    4.35521148540168,    3.86769104037540,    3.39550013091521,    0.20000000000000]
                    # chord_to_use = [2.599786083594759, 3.533583676452576, 4.44802969320108, 4.50303999525842, 4.23785803303602, 4.0690023547926, 4.355211847489540, 3.8676910403754, 3.39550013091524, 0.2]
                    chord_to_use = chord

                for k in range(self.num_blade_segments):
                    turb_i_chord.append(Constant(chord_to_use[k]))
                    turb_i_twist.append(Constant(twist[k]))
                    turb_i_lift.append(Constant(cl[k]))
                    turb_i_drag.append(Constant(cd[k]))

                self.fprint('Turbine #%d: Chord = %s' % (turb_i, np.array2string(chord_to_use, precision=6, separator=', ')))

                self.mchord.append(turb_i_chord)
                self.mtwist.append(turb_i_twist)
                self.mcl.append(turb_i_lift)
                self.mcd.append(turb_i_drag)
            self.chord = np.array(self.mchord,dtype=float)
            self.cl = np.array(self.mcl,dtype=float)
            self.cd = np.array(self.mcd,dtype=float)
            self.farm.baseline_chord = np.array(self.chord[0])/self.chord_factor

            self.cyld_expr_list = [None]*self.farm.numturbs
            self.tf_list = self.farm.CalculateActuatorLineTurbineForces(self, simTime)

            self.CopyALMtoWindFarm()
            tf = sum(self.tf_list)

        else:
            raise ValueError("Unknown turbine method: "+self.farm.turbine_method)
        
        return tf


    ############################################
    ############################################
    ### this is a hack to accommodate the fact that
    ### alm values are not stored in the wind_farm
    ### object. eventually we might want to move
    ### it there.
    ############################################
    ############################################
    def CopyALMtoWindFarm(self):
        self.farm.mcl = self.mcl
        self.farm.mcd = self.mcd
        self.farm.mchord = self.mchord
        self.farm.cl = self.cl
        self.farm.cd = self.cd
        self.farm.chord = self.chord
        self.farm.num_blade_segments = self.num_blade_segments


    def ChangeWindAngle(self,inflow_angle):
        """
        This function recomputes all necessary components for a new wind direction

        Args: 
            inflow_angle (float): The new wind angle in radians
        """
        adj_start = time.time()
        self.fprint("Adjusting Wind Angle",special="header")
        self.fprint("New Angle: {:1.8f} rads".format(inflow_angle))
        self.dom.RecomputeBoundaryMarkers(inflow_angle)
        self.bd.RecomputeVelocity(inflow_angle)
        self.ComputeFunctional(inflow_angle)

        adj_stop = time.time()
        self.fprint("Wind Angle Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")

        # self.tf.assign(tf_temp)

    def ChangeWindSpeed(self,inflow_speed):
        adj_start = time.time()
        self.fprint("Adjusting Wind Speed",special="header")
        self.fprint("New Speed: {:1.8f} m/s".format(inflow_speed/self.dom.xscale))
        self.bd.HH_vel = inflow_speed
        adj_stop = time.time()
        self.fprint("Wind Speed Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")

    def UpdateActuatorLineControls(self, c_lift = None, c_drag = None, chord = None, yaw = None, turb_i = 0):

        if c_lift is not None:
            cl = np.array(c_lift, dtype = float)
            self.cl[turb_i] = cl
            for k in range(self.num_blade_segments):
                self.mcl[turb_i][k] = Constant(cl[k])
        if c_drag is not None:
            cd = np.array(c_drag, dtype = float)
            self.cd[turb_i] = cd
            for k in range(self.num_blade_segments):
                self.mcd[turb_i][k] = Constant(cd[k])
        if chord is not None:
            chord = np.array(chord, dtype = float)
            self.chord[turb_i] = chord
            for k in range(self.num_blade_segments):
                self.mchord[turb_i][k] = Constant(chord[k])
        if yaw is not None:
            yaw = float(yaw)
            self.farm.yaw[turb_i] = yaw
            self.farm.myaw[turb_i] = Constant(yaw)
        

        self.CopyALMtoWindFarm()


class StabilizedProblem(GenericProblem):
    """
    The StabilizedProblem setup everything required for solving Navier-Stokes with 
    a stabilization term

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_conditions):
        super(StabilizedProblem, self).__init__(domain,windfarm,function_space,boundary_conditions)
        
        ### Create Functional ###
        self.ComputeFunctional()
        self.DebugOutput()


    def ComputeFunctional(self,inflow_angle=None):
        self.fprint("Setting Up Stabilized Problem",special="header")

        ### Create the test/trial/functions ###
        self.up_k = Function(self.fs.W)
        self.u_k,self.p_k = split(self.up_k)
        v,q = TestFunctions(self.fs.W)

        ### Set the x scaling ###
        Sx = self.dom.xscale

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_k.assign(self.bd.u0)

        # mem0=memory_usage()[0]
        # mem_out, self.tf = memory_usage((self.ComputeTurbineForce,(self.u_k,inflow_angle),{}),max_usage=True,retval=True,max_iterations=1)
        # self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))
        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)

        ### These constants will be moved into the params file ###
        f = Constant((0.0,)*self.dom.dim)
        f.rename("f","f")
        
        nu = self.viscosity
        vonKarman=0.41
        eps=Constant(1.0)
        eps.rename("eps","eps")

        self.fprint("Viscosity:                 {:1.2e}".format(float(self.viscosity)))
        self.fprint("Max Mixing Length:         {:1.2e}".format(float(self.lmax)))
        self.fprint("Stabilization Coefficient: {:1.2e}".format(float(eps)))

        ### Calculate the stresses and viscosities ###
        S = sqrt(2*inner(0.5*(grad(self.u_k)+grad(self.u_k).T),0.5*(grad(self.u_k)+grad(self.u_k).T)))

        ### Create l_mix based on distance to the ground ###
        if self.dom.dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,(1.+np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,self.lmax)))
        else:
            l_mix = Constant(vonKarman*self.farm.HH[0]/(1+(vonKarman*self.farm.HH[0]/self.lmax)))
        # l_mix = Expression("x[2]/8.",degree=2)
        l_mix.rename("l_mix","l_mix")
        
        ### Calculate nu_T
        self.nu_T=l_mix**2.*S
        # self.nu_T=Constant(0.0)
        self.ReyStress=self.nu_T*grad(self.u_k)
        self.vertKE= self.ReyStress[0,2]*self.u_k[0]

        ### Create the functional ###
        # if self.farm.yaw[0]**2 > 1e-4:
        #     self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
        # else :
        # self.F = inner(grad(self.u_k)*self.u_k, v)*dx + Sx*Sx*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx# + inner(self.tf,v)*dx 
        self.F = inner(grad(self.u_k)*self.u_k, v)*dx + Sx*Sx*(nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
        # self.F_sans_tf =  (1.0)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx
        # self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf*(self.u_k[0]**2+self.u_k[1]**2),v)*dx 

        ### Add in the Stabilizing term ###
        # stab = - eps*inner(grad(q), grad(self.p_k))*dx - eps*inner(grad(q), dot(grad(self.u_k), self.u_k))*dx 
        stab = - eps*inner(grad(q), grad(self.p_k))*dx - eps*inner(grad(q), dot(grad(self.u_k), self.u_k))*dx 
        # stab_sans_tf = - eps*inner(grad(q), grad(self.p_k))*dx 

        self.F += stab
        # self.F_sans_tf += stab

        if self.use_25d_model and self.dom.dim == 2 :
            if self.dom.dim == 3:
                raise ValueError("The 2.5D model requires a 2D simulation.")

            self.fprint("Using 2.5D model")
            dudx = Dx(self.u_next[0], 0)
            dvdy = Dx(self.u_next[1], 1)

            if inflow_angle is None:
                term25 = dvdy*q*dx
            else:
                term25 = (abs(sin(inflow_angle))*dudx*q + abs(cos(inflow_angle))*dvdy*q)*dx

            self.F -= term25

        self.fprint("Stabilized Problem Setup",special="footer")


class TaylorHoodProblem(GenericProblem):
    """
    The TaylorHoodProblem sets up everything required for solving Navier-Stokes 

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_conditions):
        super(TaylorHoodProblem, self).__init__(domain,windfarm,function_space,boundary_conditions)

        ### Create Functional ###
        self.ComputeFunctional()
        self.DebugOutput()

    def ComputeFunctional(self,inflow_angle=None):
        self.fprint("Setting Up Taylor-Hood Problem",special="header")

        ### These constants will be moved into the params file ###
        f = Constant((0.0,)*self.dom.dim)
        vonKarman=0.41
        eps=Constant(0.01)
        nu = self.viscosity


        self.fprint("Viscosity:         {:1.2e}".format(float(self.viscosity)))
        self.fprint("Max Mixing Length: {:1.2e}".format(float(self.lmax)))

        ### Create the test/trial/functions ###
        self.up_k = Function(self.fs.W)
        self.u_k,self.p_k = split(self.up_k)
        v,q = TestFunctions(self.fs.W)

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_k.assign(self.bd.u0)

        ### Calculate the stresses and viscosities ###
        S = sqrt(2.*inner(0.5*(grad(self.u_k)+grad(self.u_k).T),0.5*(grad(self.u_k)+grad(self.u_k).T)))

        ### Create l_mix based on distance to the ground ###
        if self.dom.dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:],(1.+np.divide(vonKarman*self.bd.depth.vector()[:],self.lmax)))
        else:
            l_mix = Constant(vonKarman*self.farm.HH[0]/(1+(vonKarman*self.farm.HH[0]/self.lmax)))

        ### Calculate nu_T
        self.nu_T=l_mix**2.*S

        ### Create the turbine force ###
        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)

        ### Create the functional ###
        self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 

        if self.use_25d_model:
            if self.dom.dim == 3:
                raise ValueError("The 2.5D model requires a 2D simulation.")

            self.fprint("Using 2.5D model")
            dudx = Dx(self.u_k[0], 0)
            dvdy = Dx(self.u_k[1], 1)

            if inflow_angle is None:
                term25 = dvdy*q*dx
            else:
                term25 = (abs(sin(inflow_angle))*dudx*q + abs(cos(inflow_angle))*dvdy*q)*dx

            self.F -= term25

        self.fprint("Taylor-Hood Problem Setup",special="footer")

# ================================================================

class IterativeSteady(GenericProblem):
    """
    The IterativeSteady sets up everything required for solving Navier-Stokes using
    the SIMPLE algorithm

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self, domain, windfarm, function_space, boundary_conditions):
        super(IterativeSteady, self).__init__(domain, windfarm, function_space, boundary_conditions)
        self.fprint("Setting Up *Iterative* Steady Problem", special="header")

        ### Create Functional ###
        self.ComputeFunctional()

    def ComputeFunctional(self,inflow_angle=None):
        # ================================================================

        # Define fluid properties
        # FIXME: These should probably be set in params.yaml input filt
        # nu = 1/10000
        nu = Constant(self.viscosity)

        # Trial functions for velocity and pressure
        u = TrialFunction(self.fs.V)
        p = TrialFunction(self.fs.Q)

        # Test functions for velocity and pressure
        v = TestFunction(self.fs.V)
        q = TestFunction(self.fs.Q)

        ### Define Velocity Functions ###
        self.u_k = Function(self.fs.V, name="u_k")
        self.u_hat = Function(self.fs.V)
        self.u_s = Function(self.fs.V)
        self.du = Function(self.fs.V)
        self.u_k_old = Function(self.fs.V)

        ### Define Pressure Functions ###
        self.p_k = Function(self.fs.Q, name="p_k")
        self.p_s = Function(self.fs.Q)
        self.dp = Function(self.fs.Q)
        self.p_k_old = Function(self.fs.Q)

        U_CN = 0.5*(u + self.u_k) 

        # Adams-Bashforth velocity
        # U_AB = 1.5*u_k - 0.5*u_k_old # Time level k+1/2
        U_AB = 2.0*self.u_k - 1.0*self.u_k_old # Time level k+1

        ### Define Eddy viscosity ###
        use_eddy_viscosity = True

        if use_eddy_viscosity:
            # Define filter scale (all models use this)
            filter_scale = CellVolume(self.dom.mesh)**(1.0/self.dom.dim)

            # Strain rate tensor, 0.5*(du_i/dx_j + du_j/dx_i)
            Sij = sym(grad(U_AB))

            # sqrt(Sij*Sij)
            strainMag = (2.0*inner(Sij, Sij))**0.5

            # Smagorinsky constant, typically close to 0.17
            Cs = 0.17

            vonKarman=0.41
            # lmax Default = 15.0
            # l_mix = Constant(vonKarman*self.farm.HH[0]/(1+(vonKarman*self.farm.HH[0]/self.lmax)))

            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:],(1.+np.divide(vonKarman*self.bd.depth.vector()[:],self.lmax)))

            # Eddy viscosity
            # self.nu_T = Cs**2 * filter_scale**2 * strainMag
            self.nu_T = l_mix**2 * strainMag

        else:
            self.nu_T = Constant(0.0)

        # ================================================================

        # FIXME: This up_k function is only present to avoid errors  
        # during assignments in GenericSolver.__init__

        # Create the combined function space
        self.up_k = Function(self.fs.W)

        # Create the turbine force
        # FIXME: Should this be set by a numpy array operation or a fenics function?
        # self.tf = self.farm.TurbineForce(self.fs, self.dom.mesh, self.u_k2)
        # self.tf = Function(self.fs.V)

        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)
        # self.u_k.assign(self.bd.bc_velocity)

        # self.u_k2.vector()[:] = 0.0
        # self.u_k1.vector()[:] = 0.0

        # ================================================================

        # Solve for u_hat, a velocity estimate which doesn't include pressure gradient effects
        F1 = inner(dot(self.u_k, nabla_grad(u)), v)*dx \
           + (nu+self.nu_T)*inner(grad(u), grad(v))*dx \
           + inner(self.tf, v)*dx 

        self.F1_lhs = lhs(F1)
        self.F1_rhs = rhs(F1)


        # Use u_hat to solve for the pressure field
        self.F2_lhs = inner(grad(p), grad(q))*dx
        self.F2_rhs = - div(self.u_hat)*q*dx

        self.dt_1 = Constant(1) # Becomes (1/dt, 0) for (unsteady, steady) state
        self.dt_2 = Constant(1) # Becomes (1/dt, 1) for (unsteady, steady) state
        self.dt_3 = Constant(1) # Becomes (  dt, 1) for (unsteady, steady) state

        # Solve for u_star, a predicted velocity which includes the pressure gradient
        F3 = inner(dot(self.u_k, nabla_grad(u)), v)*dx \
           + (nu+self.nu_T)*inner(grad(u), grad(v))*dx \
           + inner(self.tf, v)*dx \
           + inner(grad(self.p_k), v)*dx \
           + self.dt_1*inner(u - self.u_k, v)*dx


        self.F3_lhs = lhs(F3)
        self.F3_rhs = rhs(F3)


        # Define variational problem for step 2: pressure correction
        self.F4_lhs = inner(grad(p), grad(q))*dx
        self.F4_rhs = - self.dt_2*div(self.u_s)*q*dx


        # Define variational problem for step 3: velocity update
        self.F5_lhs = inner(u, v)*dx
        self.F5_rhs = - self.dt_3*inner(grad(self.dp), v)*dx

        # ================================================================

        self.fprint("*Iterative* Steady Problem Setup",special="footer")

# ================================================================

class UnsteadyProblem(GenericProblem):
    """
    The UnsteadyProblem sets up everything required for solving Navier-Stokes using
    a fractional-step method with an adaptive timestep size

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self, domain, windfarm, function_space, boundary_conditions):
        super(UnsteadyProblem, self).__init__(domain, windfarm, function_space, boundary_conditions)
        self.fprint("Setting Up Unsteady Problem", special="header")

        ### Create Functional ###
        self.ComputeFunctional()
        self.DebugOutput()

    def ComputeFunctional(self,inflow_angle=None):
        # ================================================================

        # Define fluid properties
        # FIXME: These should probably be set in params.yaml input filt
        # nu = 1/10000
        rho = 1
        nu_c = Constant(self.viscosity)
        rho_c = Constant(rho)

        # Define time step size (this value is used only for step 1 if adaptive timestepping is used)
        # FIXME: change variable name to avoid confusion within dolfin adjoint
        self.dt = 0.1*self.dom.mesh.hmin()/self.bd.HH_vel
        self.dt = 0.05
        self.dt_c  = Constant(self.dt)

        self.fprint("Viscosity: {:1.2e}".format(float(self.viscosity)))
        self.fprint("Density:   {:1.2e}".format(float(rho)))

        # Define trial and test functions for velocity
        u = TrialFunction(self.fs.V)
        v = TestFunction(self.fs.V)

        # Define trial and test functions for pressure
        p = TrialFunction(self.fs.Q)
        q = TestFunction(self.fs.Q)

        # Define functions for velocity solutions
        # >> _k = current (step k)
        # >> _k1 = previous (step k-1)
        # >> _k2 = double previous (step k-2)
        self.u_k = Function(self.fs.V)
        self.u_k1 = Function(self.fs.V)
        self.u_k2 = Function(self.fs.V)

        # Seed previous velocity fields with the chosen initial condition
        self.u_k.assign(self.bd.bc_velocity)
        self.u_k1.assign(self.bd.bc_velocity)
        self.u_k2.assign(self.bd.bc_velocity)

        # Calculate Reynolds stress 
        self.uk_sum = Function(self.fs.V)
        self.uk_sum.assign(self.dt_c*self.u_k)
        self.vertKE = Function(self.fs.Q)

        # Define functions for pressure solutions
        # >> _k = current (step k)
        # >> _k1 = previous (step k-1)
        self.p_k  = Function(self.fs.Q)
        self.p_k1 = Function(self.fs.Q)

        # Seed previous pressure fields with the chosen initial condition
        self.p_k1.assign(self.bd.bc_pressure)

        # ================================================================

        # Crank-Nicolson velocity
        U_CN  = 0.5*(u + self.u_k1)

        # Adams-Bashforth projected velocity
        U_AB = 1.5*self.u_k1 - 0.5*self.u_k2

        # ================================================================

        # Calculate eddy viscosity, if using
        use_eddy_viscosity = True

        if use_eddy_viscosity:
            # Define filter scale
            filter_scale = CellVolume(self.dom.mesh)**(1.0/self.dom.dim)

            # Strain rate tensor, 0.5*(du_i/dx_j + du_j/dx_i)
            Sij = sym(nabla_grad(U_AB))

            # sqrt(Sij*Sij)
            strainMag = (2.0*inner(Sij, Sij))**0.5

            # Smagorinsky constant, typically around 0.17
            Cs = 0.17

            # Eddy viscosity
            self.nu_T = Cs**2 * filter_scale**2 * strainMag
        else:
            self.nu_T = Constant(0)

        # ================================================================

        # FIXME: This up_k function is only present to avoid errors  
        # during assignments in GenericSolver.__init__

        # Create the combined function space
        self.up_k = Function(self.fs.W)

        # Create the turbine force
        # FIXME: Should this be set by a numpy array operation or a fenics function?
        # self.tf = self.farm.TurbineForce(self.fs, self.dom.mesh, self.u_k2)
        # self.tf = Function(self.fs.V)

        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)
        self.u_k.assign(self.bd.bc_velocity)


        # self.u_k2.vector()[:] = 0.0
        # self.u_k1.vector()[:] = 0.0

        # ================================================================

        # Define variational problem for step 1: tentative velocity
        # F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
        #    + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
        #    + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
        #    + dot(nabla_grad(self.p_k1), v)*dx \
        #    - dot(-self.tf, v)*dx

        # F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
        #    + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
        #    + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
        #    + dot(nabla_grad(self.p_k1), v)*dx \
        #    - dot(self.tf, v)*dx

        F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
           + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
           + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
           + inner(grad(self.p_k1), v)*dx \
           - dot(self.tf, v)*dx

        self.a1 = lhs(F1)
        self.L1 = rhs(F1)

        # Define variational problem for step 2: pressure correction
        # self.a2 = dot(nabla_grad(p), nabla_grad(q))*dx
        # self.L2 = dot(nabla_grad(self.p_k1), nabla_grad(q))*dx - (1.0/self.dt_c)*div(self.u_k)*q*dx
        self.a2 = inner(grad(p), grad(q))*dx
        self.L2 = inner(grad(self.p_k1), grad(q))*dx - (1.0/self.dt_c)*div(self.u_k)*q*dx

        # phi = p - self.p_k
        # F2 = inner(grad(q), grad(phi))*dx - (1.0/self.dt_c)*div(u_k)*q*dx
        # self.a2 = lhs(F2)
        # self.L2 = rhs(F2)

        # Define variational problem for step 3: velocity update
        # self.a3 = dot(u, v)*dx
        # self.L3 = dot(self.u_k, v)*dx - self.dt_c*dot(nabla_grad(self.p_k - self.p_k1), v)*dx
        self.a3 = inner(u, v)*dx
        self.L3 = inner(self.u_k, v)*dx - self.dt_c*inner(grad(self.p_k - self.p_k1), v)*dx

        # F3 = inner(u, v)*dx - inner(self.u_k, v)*dx + self.dt_c*inner(phi, v)*dx
        # self.a3 = lhs(F3)
        # self.L3 = rhs(F3)


        # ================================================================

        self.fprint("Unsteady Problem Setup",special="footer")

# # ================================================================

#     def UpdateActuatorLineControls(self, c_lift = None, c_drag = None):

#         if c_lift is not None:
#             cl = np.array(c_lift, dtype = float)
#         if c_drag is not None:
#             cd = np.array(c_drag, dtype = float)

#         for k in range(self.num_blade_segments):
#             self.mcl[k] = Constant(cl[k])
#             self.mcd[k] = Constant(cd[k])

# # ================================================================

