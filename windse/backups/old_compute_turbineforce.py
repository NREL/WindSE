OLD COMPUTE TURBINE FORCE





    def ComputeTurbineForce(self,u,inflow_angle,simTime=0.0):

        ### Compute the relative yaw angle ###
        # if inflow_angle is not None:
        #     inflow_angle = inflow_angle-self.dom.inflow_angle
        # else:
        #     inflow_angle = self.dom.inflow_angle
        if inflow_angle is None:
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
            # if not os.path.exists(aoa_folder): os.makedirs(aoa_folder)
            # if not os.path.exists(force_folder): os.makedirs(force_folder)
            if not os.path.exists(aoa_folder) and self.params.rank==0: os.makedirs(aoa_folder)
            if not os.path.exists(force_folder) and self.params.rank==0: os.makedirs(force_folder)
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
