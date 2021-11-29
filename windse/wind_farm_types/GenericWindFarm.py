from windse import windse_parameters

class GenericWindFarm(object):
    """
    A GenericWindFarm contains on the basic functions and attributes required by all wind farm objects.
    
    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self, dom):
        """
        Store anything needed prior to setup
        """

        # Store WindSE Objects
        self.dom = dom
        self.params = windse_parameters
        self.fprint = self.params.fprint
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode

        # Load any global parameter
        self.turbine_type = self.params["turbines"]["type"]

        # Init blank turbine list
        self.turbines = []

        # Setup the wind farm (eventually this will happen outside the init)
        self.fprint(f"Generating {self.name}",special="header")
        self.Setup()
        self.fprint(f"{self.name} Generated",special="footer")

    def Setup(self):
        """
        This function builds the wind farm as well as sets up the turbines
        """  
        self.LoadParameters()      
        self.initial_turbine_locations = self.InitializeTurbineLocations()

        self.fprint("Setting Up Turbines",special="header")
        self.SetupTurbines()
        self.fprint("Turbines Set up",special="footer")

        self.DebugOutput() 

    def LoadParameters(self):
        """
        This function will parse the parameters from the yaml file
        """  
        raise NotImplementedError(type(self))

    def InitializeTurbineLocations():
        """
        This function will compute the initial locations of the turbines in the farm
        It must return an Nx2 array where N is the number of turbines
        """
        raise NotImplementedError(type(self))

    def SetupTurbines(self):
        """
        Using the parameters and initial locations this function will populate the list of turbines
        """
        from windse.turbine_types import turbine_dict
        turbine_method = turbine_dict[self.turbine_type]
        for i,(x,y) in enumerate(self.initial_turbine_locations):
            self.turbines.append(turbine_method(i,x,y,self.dom))

    def ComputeTurbineForce(self):
        """
        Iterates over the turbines and adds up the turbine forces
        """
        pass

    def UpdateControls(self):
        """
        """
        pass

    def UpdateTurbines(self):
        """
        Updates the turbines
        """
        pass

    def UpdateTurbineForce(self):
        """
        Updates the turbine force
        """
        pass

    def GetTurbineLocation(self):
        """
        returns a nx3 numpy array containing the x, y, and z location of the turbine hubs
        """
        pass

    def DebugOutput(self):
        """
        This function computes and save any output needed for regression tests
        """
        pass

    def PlotFarm(self):
        """
        saves a plot of all the turbines current location and yaw
        """
        pass

    def PlotChord(self):
        """
        saves a plot of the chord for each turbine?
        """
        pass

    def SaveWindFarm(self):
        """
        saves a text file containing the wind_farm and turbine parameters
        """
        pass

    def SaveTurbineForce(self):
        """
        save the turbine force as a paraviewable function
        """
        pass




















    ########################################################################################################
    ############## It would be smart to eventually move these to a separate refinement module ##############
    ########################################################################################################
 
    def SimpleRefine(self,radius,expand_factor=1):
        if self.numturbs == 0:
            return

        self.fprint("Cylinder Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)
        
        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_x = np.array(self.x)
        turb_y = np.array(self.y)
        if self.dom.dim == 3:
            turb_z0 = self.dom.z_range[0]-radius
            turb_z1 = self.z+radius

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Find the minimum distance for each turbine with the vertices ###
            x_diff = np.power(np.subtract.outer(x,turb_x),2.0)
            y_diff = np.power(np.subtract.outer(y,turb_y),2.0)
            min_r = np.min(x_diff+y_diff,axis=0)

            ### Determine if cell is in radius for each turbine ###
            in_circle = min_r <= radius**2.0

            ### Determine if in cylinder ###
            if d == 3:
                in_z = np.logical_and(turb_z0 <= max(z), turb_z1 >= min(z))
                near_turbine = np.logical_and(in_circle, in_z)
            else:
                near_turbine = in_circle

            ### mark if cell is near any cylinder ###
            if any(near_turbine):
                cell_f[cell] = True

        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.CalculateHeights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")


    def WakeRefine(self,radius,length,theta=0.0,expand_factor=1,centered=False):
        if self.numturbs == 0:
            return

        self.fprint("Wake Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius/2.0
        length = length+2*(expand_factor-1)*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)

        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_x = np.array(self.x)
        turb_y = np.array(self.y)
        if self.dom.dim == 3:
            turb_z = np.array(self.z)

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Rotate the Cylinder about the turbine axis ###
            x_diff = np.subtract.outer(x,turb_x)
            y_diff = np.subtract.outer(y,turb_y)
            x = (np.cos(theta)*(x_diff)-np.sin(theta)*(y_diff) + turb_x)
            y = (np.sin(theta)*(x_diff)+np.cos(theta)*(y_diff) + turb_y)

            ### Determine if in wake ###
            if centered:
                # Center the refinement region around the turbine
                # upstream and downstream by length/2
                x0 = turb_x - length/2.0
                x1 = turb_x + length/2.0
            else:
                # Otherwise, refine the default amount upstream (R)
                # and the full length amount downstreeam
                x0 = turb_x - radius
                x1 = turb_x + length

            in_wake = np.logical_and(np.greater(x,x0),np.less(x,x1))
            in_wake = np.any(in_wake,axis=0)

            ### Find the minimum distance for each turbine with the vertices ###
            y_diff = y-turb_y
            if d == 3:
                z_diff = np.subtract.outer(z,turb_z)
                min_r = np.min(np.power(y_diff,2.0)+np.power(z_diff,2.0),axis=0)
            else:
                min_r = np.min(np.power(y_diff,2.0),axis=0)

            ### Determine if cell is in radius for each turbine ###
            in_circle = min_r <= radius**2.0
            near_turbine = np.logical_and(in_circle, in_wake)

            ### mark if cell is near any cylinder ###
            if any(near_turbine):
                cell_f[cell] = True


        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.CalculateHeights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def TearRefine(self,radius,theta=0.0,expand_factor=1):
        if self.numturbs == 0:
            return

        self.fprint("Tear Drop Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)
        
        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_x = np.array(self.x)
        turb_y = np.array(self.y)
        if self.dom.dim == 3:
            turb_z0 = self.dom.z_range[0]-radius
            turb_z1 = self.z+radius

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Rotate the Cylinder about the turbine axis ###
            x_diff = np.subtract.outer(x,turb_x)
            y_diff = np.subtract.outer(y,turb_y)
            x = (np.cos(theta)*(x_diff)-np.sin(theta)*(y_diff) + turb_x)
            y = (np.sin(theta)*(x_diff)+np.cos(theta)*(y_diff) + turb_y)
            x_diff = x - turb_x
            y_diff = y - turb_y

            ### Find Closest Turbine ###
            min_dist = np.min(np.power(x_diff,2.0)+np.power(y_diff,2.0),axis=0)
            min_turb_id = np.argmin(min_dist)

            # ### Do something based on upstream or downstream ###
            # if min(x[:,min_turb_id]-turb_x[min_turb_id]) <= 0:
            #     val = -min_turb_id-1
            # else:
            #     val = min_turb_id+1

            # ### Determine if in z_range ###
            # if d == 3:
            #     in_z = turb_z0 <= max(z) and turb_z1[min_turb_id] >= min(z)
            #     if in_z:
            #         near_turbine = val
            #     else:
            #         near_turbine = 0
            # else:
            #     near_turbine = val

            ### Check if upstream or downstream and adjust accordingly ###
            cl_x = turb_x[min_turb_id]
            cl_y = turb_y[min_turb_id]
            dx = min(x[:,min_turb_id]-cl_x)
            dy = min(y[:,min_turb_id]-cl_y)
            if dx <= 0:
                dist = (dx*1.5)**2 + dy**2
            else:
                dist = (dx/2)**2 + dy**2

            ### determine if in z-range ###
            if d == 3:
                in_z = turb_z0 <= max(z) and turb_z1[min_turb_id] >= min(z)
                near_turbine = in_z and dist <= radius**2
            else:
                near_turbine = dist <= radius**2
            ### mark if cell is near any cylinder ###
            if near_turbine:
                cell_f[cell] = True


        # File("test.pvd") << cell_f
        # exit()

        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.CalculateHeights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

    def SphereRefine(self,radius,expand_factor=1):
        if self.numturbs == 0:
            return

        self.fprint("Sphere Refinement Near Turbines",special="header")
        refine_start = time.time()

        ### Calculate expanded values ###
        radius = expand_factor*radius

        ### Create the cell markers ###
        cell_f = MeshFunction('bool', self.dom.mesh, self.dom.mesh.geometry().dim(),False)
        
        ### Get Dimension ###
        n = self.numturbs
        d = self.dom.dim

        ### Get Turbine Coordinates ###
        turb_x = np.array(self.x)
        turb_y = np.array(self.y)
        if self.dom.dim == 3:
            turb_z = np.array(self.z)
            

        self.fprint("Marking Near Turbine")
        mark_start = time.time()
        for cell in cells(self.dom.mesh):
            ### Get Points of all cell vertices ##
            pt = cell.get_vertex_coordinates()
            x = pt[0::d]
            y = pt[1::d]
            if d == 3:
                z = pt[2::d]

            ### Find the minimum distance for each turbine with the vertices ###
            min_r  = np.power(np.subtract.outer(x,turb_x),2.0)
            min_r += np.power(np.subtract.outer(y,turb_y),2.0)
            if d == 3:
                min_r += np.power(np.subtract.outer(z,turb_z),2.0)
            min_r = np.min(min_r,axis=0)

            ### Determine if cell is in radius for each turbine ###
            in_sphere = min_r <= radius**2.0

            ### mark if cell is near any cylinder ###
            if any(in_sphere):
                cell_f[cell] = True

        mark_stop = time.time()
        self.fprint("Marking Finished: {:1.2f} s".format(mark_stop-mark_start))

        ### Refine Mesh ###
        self.dom.Refine(cell_f)

        ### Recompute Heights with new mesh ###
        self.CalculateHeights()

        refine_stop = time.time()
        self.fprint("Mesh Refinement Finished: {:1.2f} s".format(refine_stop-refine_start),special="footer")

        