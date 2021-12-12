# import windse data
from windse import windse_parameters
if windse_parameters.dolfin_adjoint:
    from windse.blocks import blockify, BaseHeightBlock

# import dolfin functions
from . import cos, sin, Constant

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

        # blockify custom functions so dolfin adjoint can track them
        if self.params.performing_opt_calc:
            self.base_height = blockify(self.base_height,BaseHeightBlock,block_kwargs={"ground":self.dom.Ground})

        # complete the setup
        self.setup()

    def setup(self):
        """
        This function takes the init data and set it up
        """  
        self.load_parameters()
        self.create_controls()
        self.calculate_heights()
        self.debug_output() 

        self.fprint("Turbine Hub {:d}: ({: 1.5f}, {: 1.5f}, {: 1.5f})".format(self.index, self.x, self.y, self.z))

    def load_parameters(self):
        """
        This function will parse the parameters from the yaml file
        """  
        raise NotImplementedError(type(self))

    def create_controls(self):
        """
        convert any initial data to dolfin constants of the form self.m<name>, where self.<name> is already defined
        additionally define the self.controls_list which is just a list of strings [<name1>,<name1>] identical to the forms above

        """
        raise NotImplementedError(type(self))

    def update_controls(self):
        """
        updates the controls during optimzation so that the output functions display the correct information
        """
        raise NotImplementedError(type(self))

    def base_height(self,x,y,dx=0,dy=0):
        """
        calculate the high of the ground at (x,y) based on the domains Ground function
        """
        val = self.dom.Ground(float(x),float(y),dx=dx,dy=dy)
        return Constant(val)

    def calculate_heights(self):
        """
        This function calculates the absolute heights of each turbine.
        """

        # use the ground function from the domain co calculate the height of the turbine hub (this preserves dolfin_adjoint)
        self.mz = self.base_height(self.mx,self.my)+float(self.HH)

        # store the float for easy access
        self.z = float(self.mz)

        # store the height of the ground 
        self.ground = self.z - self.HH

    def yaw_turbine(self,x,x0,yaw):
        """
        This function yaws the turbines when creating the turbine force.

        Args:
            x (dolfin.SpatialCoordinate): the space variable, x
            x0 (list): the location of the turbine to be yawed
            yaw (float): the yaw value in radians
        """
        xrot =   cos(yaw)*(x[0]-x0[0]) + sin(yaw)*(x[1]-x0[1])
        yrot = - sin(yaw)*(x[0]-x0[0]) + cos(yaw)*(x[1]-x0[1])
        if self.dom.dim == 3:
            zrot = x[2]-x0[2]
        else:
            zrot = 0.0
        
        return [xrot,yrot,zrot]

    def rotate_tower(self):
        """
        rotates the turbine by rotating the x,y,z coordinates of the domain
        """
        pass

    def translate_tower(self):
        """
        translates the turbine by shifting the x,y,z coordinates of the domain
        """
        pass

    def force(self):
        """
        computes the turbine force
        """
        raise NotImplementedError(type(self))

    def force_gradient(self):
        """
        will give the gradient of the force with respect to a given control
        """
        raise NotImplementedError(type(self))

    def power(self):
        """
        computes the power produced by this turbine
        """
        raise NotImplementedError(type(self))

    def power_gradient(self):
        """
        will give the gradient of the power with respect to a given control
        """
        raise NotImplementedError(type(self))

    def debug_output(self):
        """
        This function computes and save any output needed for regression tests
        """
        pass

    def prepare_saved_functions(self, func_list):
        '''
        This method prepares any dolfin function that the specific turbine type would like to save. 
        func_list is the output file, but it is modified in place. it has the format [[func1,"func1_name"],[func2,"func2_name"],...]
        '''
        raise NotImplementedError(type(self))

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