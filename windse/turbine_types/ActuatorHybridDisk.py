# import windse function
from . import GenericTurbine

# import dolfin functions
from . import Constant, SpatialCoordinate, as_vector, cos, sin, exp, sqrt, dot, project, assemble, dx


# import other modules
import numpy as np
from scipy.special import gamma
import scipy.interpolate as interp

class ActuatorHybridDisk(GenericTurbine):

    def __init__(self, i,x,y,dom,imported_params=None):
        # Define the acceptable column of an wind_farm.csv imput file        
        self.yaml_inputs = ["HH", "RD", "yaw", "thickness", "axial", "force"]

        # Init turbine
        super(ActuatorHybridDisk, self).__init__(i,x,y,dom,imported_params)
    

    def load_parameters(self):
        self.HH             = self.params["turbines"]["HH"]
        self.RD             = self.params["turbines"]["RD"]
        self.yaw            = self.params["turbines"]["yaw"]
        self.thickness      = self.params["turbines"]["thickness"]
        self.axial          = self.params["turbines"]["axial"]
        self.force          = self.params["turbines"]["force"]
        self.blade_segments = self.params["turbines"]["blade_segments"]
        self.read_turb_data = self.params["turbines"]["read_turb_data"]
    
    def create_controls(self, initial_call_from_setup=True):

        if initial_call_from_setup:
            self.mx     = Constant(self.x, name="x_{:d}".format(self.index))
            self.my     = Constant(self.y, name="y_{:d}".format(self.index))
            self.myaw   = Constant(self.yaw, name="yaw_{:d}".format(self.index))
            self.maxial = Constant(self.axial, name="axial_{:d}".format(self.index))
        else:
            self.mthrust = []
            self.mspin = []
            for i in range(self.blade_segments):
                self.mthrust.append(Constant(self.thrust[i], name=f"thrust_({self.index},{i})"))
                self.mspin.append(Constant(self.spin[i], name=f"spin_({self.index},{i})"))

        # The control list is very important for extracting information from the optimzation step.
        self.controls_list = ["x","y","yaw","axial","thrust","spin"] 

    def init_blade_properties(self):

        # Aliasing
        if isinstance(self.blade_segments,str):  
            self.fprint(f'Blade segments set to "{self.blade_segments}," using 4 segments by default')
            self.blade_segments = 4
        self.num_actuator_nodes = self.blade_segments

        # Import blade data
        if self.read_turb_data is None:
            self.fprint('No turbine profile provide, using defaults')
            # These values are located at the center of the ring
            self.thrust = np.ones(self.blade_segments)
            self.spin   = np.zeros(self.blade_segments)
            # self.thrust = [1.000, 1.000, 1.000, 1.000]
            # self.spin   = [0.000, 0.000, 0.000, 0.000]
            # self.thrust = [0.600, 0.700, 1.240, 1.050]
            # self.spin   = [0.01, 0.05, 0.10, 0.07]
            # self.thrust = [0.43824640652668556, 0.5909172230419787, 0.6308489720927284, 0.6648589641278763]
            # self.spin   = [0.061621858574383126, 0.1588757152218311, 0.15187482621601978, 0.1256821489866074]
        else: 
            self.fprint('Setting thrust and spin from file \'%s\'' % (self.read_turb_data))
            actual_turbine_data = np.genfromtxt(self.read_turb_data, delimiter = ',', skip_header = 1)
            actual_x      = actual_turbine_data[:, 0]
            actual_thrust = actual_turbine_data[:, 1]
            actual_spin   = actual_turbine_data[:, 2]

            # Create interpolators for chord, lift, and drag
            interp_thrust = interp.interp1d(actual_x, actual_thrust)
            interp_spin   = interp.interp1d(actual_x, actual_spin)

            # Construct the points at which to generate interpolated values
            interp_points = np.linspace(0.0, 1.0, self.num_actuator_nodes+1)

            # Shift by h/2
            interp_points = interp_points[0:-1]+1/(2*self.num_actuator_nodes)

            # Generate the interpolated values
            self.thrust = interp_thrust(interp_points)
            self.spin   = interp_spin(interp_points)
            self.fprint(f'Blade locations: {interp_points}',offset=1)
            self.fprint(f'Thrust Values:   {self.thrust}',offset=1)
            self.fprint(f'Spin Values:     {self.spin}',offset=1)

        self.create_controls(initial_call_from_setup=False)


    def build_actuator_disk(self,inflow_angle):

        ### Alias useful values ###
        x0 = [self.mx,self.my,self.mz]
        yaw = self.myaw+self.inflow_angle
        W = self.thickness*1.0
        R = self.RD/2.0
        ma = self.maxial
        C_tprime = 4*ma/(1-ma)
        x=SpatialCoordinate(self.dom.mesh)

        ### Set up some constants ###
        S_norm = (2.0+np.pi)/(2.0*np.pi)
        T_norm = 2.0*gamma(7.0/6.0)
        A = np.pi*R**2.0 
        D_norm = np.pi*gamma(4.0/3.0)

        ### Rotate and Shift the Turbine ###
        xs = self.yaw_turbine(x,x0,yaw)

        ### Create the function that represents the Thickness of the turbine ###
        T = exp(-pow((xs[0]/W),6.0))

        ### Create the function that represents the Disk of the turbine
        r = sqrt(xs[1]**2.0+xs[2]**2.0)

        ### Calculate normalization constant ###
        volNormalization = T_norm*D_norm*W*R**(self.dom.dim-1)

        ### create range of radii
        Rs = np.linspace(0,R,self.num_actuator_nodes+1)

        # Build disk
        actuator_disk = 0
        for i in range(self.num_actuator_nodes): 
            # Build Ring
            if i==0:
                D = exp(-pow(r/Rs[i+1],6.0))
            else:
                D = exp(-pow(r/Rs[i+1],6.0))-exp(-pow(r/Rs[i],6.0))

            ### Create the function that represents the force ###
            F = -0.5*A*C_tprime*self.mthrust[i]

            ### Create the rotational force function
            RF_y = -self.mspin[i]*xs[2]/(r)
            RF_z =  self.mspin[i]*xs[1]/(r)

            ### Create direction vector
            WTGbase = as_vector((cos(yaw)-RF_y*sin(yaw),sin(yaw)+RF_y*cos(yaw),RF_z))            
            # WTGbase = as_vector((F*cos(yaw)-RF_y*sin(yaw),F*sin(yaw)+RF_y*cos(yaw),RF_z))  # this feel more right but completely fails          

            ### assemble disk
            actuator_disk += F*T*D*WTGbase/volNormalization
            # actuator_disk += T*D*WTGbase/volNormalization








        # ### Build polynomials ###
        # thrust = 0.0
        # spin = 0.0
        # for i in range(self.num_actuator_nodes):
        #     thrust += self.mthrust[i]*r**i
        #     spin  += self.mspin[i]*r**i
        # # thrust *= -(r-0.0)*(r-1.0)
        # # spin  *= -(r+0.1)*(r-1.1)

        # ### Create the function that represents the force ###
        # F = -0.5*A*C_tprime*thrust

        # 1,-z,y


        # ### Create the rotational force function
        # RF_y = -spin*xs[2]/(R*r)
        # RF_z =  spin*xs[1]/(R*r)

        # ### Create direction vector
        # WTGbase = as_vector((cos(yaw)-RF_y*sin(yaw),sin(yaw)+RF_y*cos(yaw),RF_z))
            


        # # compute disk averaged velocity in yawed case and don't project
        # actuator_disks=F_i*T*D_i*WTGbase_i/volNormalization

        return actuator_disk

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

        ### Store function space ###
        self.fs = fs
        self.inflow_angle.assign(inflow_angle)

        ### Initialize some blade parameters ###
        self.init_blade_properties()

        ### compute the space kernal and radial force
        self.actuator_disk = self.build_actuator_disk(inflow_angle)

        ### Expand the dot product
        yaw = self.myaw+self.inflow_angle
        self.tf1 = self.actuator_disk * cos(yaw)**2
        self.tf2 = self.actuator_disk * sin(yaw)**2
        self.tf3 = self.actuator_disk * 2.0 * cos(yaw) * sin(yaw)

        ### Compose full turbine force
        self.tf = self.tf1*u[0]**2+self.tf2*u[1]**2+self.tf3*u[0]*u[1]

        return self.tf

    def force_gradient(self):
        pass

    def power(self, u, inflow_angle):
        return dot(-self.tf,u)/1.0e6

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









    #### These two chord based functions are still very experimental 
    def radial_chord_force(r,chord):

        cx = np.linspace(-0.1,1.25,len(chord))
        # cx = np.insert(cx,0,-0.1)
        # cx = np.insert(cx,len(cx),1.1)
        # chord = [0.0] + chord + [0.0]

        # chord = np.array(chord,dtype=float)
        # print(chord)
        # r = np.linspace(-0.1,1.25,1000)

        force = 0
        int_force = 0
        for i in range(len(chord)):
            Li = 1.0 
            for j in range(len(chord)):
                if i!=j:
                    Li *= (r - cx[j])/(cx[i]-cx[j])

            force += chord[i]*Li
            if i !=0:
                int_force += (cx[i]-cx[i-1])*(chord[i]+chord[i-1])/2.0

        # plt.plot(r,force)
        # plt.scatter(cx,chord)
        # plt.show()
        # exit()

        return force/int_force

    def _setup_chord_force(self):
        ### this section of code is a hack to get "chord"-type disk representation ###
        if self.mchord is not None:
            if self.chord is not None:
                if self.blade_segments == "computed":
                    self.num_actuator_nodes = 10 ##### FIX THIS ####
                    self.blade_segments = self.num_actuator_nodes
                else:
                    self.num_actuator_nodes = self.blade_segments

                if self.read_turb_data:
                    print('Num blade segments: ', self.num_actuator_nodes)
                    turb_data = self.params["wind_farm"]["read_turb_data"]
                    self.fprint('Setting chord from file \'%s\'' % (turb_data))
                    actual_turbine_data = np.genfromtxt(turb_data, delimiter = ',', skip_header = 1)
                    actual_x = actual_turbine_data[:, 0]
                    actual_chord = self.chord_factor*actual_turbine_data[:, 1]
                    chord_interp = interp.interp1d(actual_x, actual_chord)
                    interp_points = np.linspace(0.0, 1.0, self.blade_segments)
                    # Generate the interpolated values
                    self.chord = chord_interp(interp_points)
                else:
                    self.chord = np.ones(self.blade_segments)
            self.num_actuator_nodes = self.blade_segments
            self.baseline_chord = copy.copy(self.chord)

            self.cl = np.ones(self.blade_segments)
            self.cd = np.ones(self.blade_segments)
            self.mcl = []
            self.mcd = []
            self.mchord = []
            for i in range(self.numturbs):
                self.mcl.append([])
                self.mcd.append([])
                self.mchord.append([])
                for j in range(self.blade_segments):
                    self.mcl[i].append(Constant(self.cl[j]))
                    self.mcd[i].append(Constant(self.cd[j]))
                    self.mchord[i].append(Constant(self.chord[j]))
            self.cl = np.array(self.mcl,dtype=float)
            self.cd = np.array(self.mcd,dtype=float)
            self.chord = np.array(self.mchord,dtype=float)

