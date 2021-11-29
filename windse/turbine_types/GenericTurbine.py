from windse import windse_parameters, BaseHeight

class GenericTurbine(object):
    """
    A GenericTurbine contains on the basic functions and attributes required by all turbine objects.
    
    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self, i,x,y,dom):
        """
        Store anything needed prior to setup
        """

        # Store windse objects
        self.dom = dom
        self.params = windse_parameters
        self.fprint = self.params.fprint
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode

        # Store turbine properties
        self.index = i
        self.x = x
        self.y = y

        self.fprint("Turbine {:d}: ({:1.5f}, {:1.5f})".format(self.index, self.x, self.y))

        self.Setup()

    def Setup(self):
        """
        This function takes the init data and set it up
        """  
        self.LoadParameters()
        self.CreateControls()
        self.CalculateHeights()
        self.DebugOutput() 

    def LoadParameters():
        """
        This function will parse the parameters from the yaml file
        """  
        raise NotImplementedError(type(self))

    def CreateControls():
        """
        convert any initial data to dolfin constants of the form self.m<name>
        """
        raise NotImplementedError(type(self))

    def UpdateControls():
        """
        updates the controls during optimzation so that the output functions display the correct information
        """
        raise NotImplementedError(type(self))

    def CalculateHeights(self):
        """
        This function calculates the absolute heights of each turbine.
        """

        # use the ground function from the domain co calculate the height of the turbine hub (this preserves dolfin_adjoint)
        self.mz = BaseHeight(self.mx,self.my,self.dom.Ground)+float(self.HH)

        # store the float for easy access
        self.z = float(self.mz)

        # store the height of the ground 
        self.ground = self.z - self.HH


    def Yaw():
        """
        yaws the rotor by rotating the x,y,z coordinates of the domain
        """
        pass

    def Rotate():
        """
        rotates the turbine by rotating the x,y,z coordinates of the domain
        """
        pass

    def Translate():
        """
        translates the turbine by shifting the x,y,z coordinates of the domain
        """
        pass

    def Force():
        """
        computes the turbine force
        """
        raise NotImplementedError(type(self))

    def ForceGradient():
        """
        will give the gradient of the force with respect to a given control
        """
        raise NotImplementedError(type(self))

    def Power():
        """
        computes the power produced by this turbine
        """
        raise NotImplementedError(type(self))

    def PowerGradient():
        """
        will give the gradient of the power with respect to a given control
        """
        raise NotImplementedError(type(self))

    def DebugOutput(self):
        """
        This function computes and save any output needed for regression tests
        """
        pass


# class TurbineForceBlock(Block):
#     """
#     A pyadjoint block that help dolfin-adjoint compute derivatives with respect to the turbine force
#     """
#     pass

# class TurbinePowerBlock(Block):
#     """
#     A pyadjoint block that help dolfin-adjoint compute derivatives with respect to the turbine force
#     """
#     pass