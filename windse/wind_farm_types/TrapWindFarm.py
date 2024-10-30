from . import GenericWindFarm
import numpy as np

class TrapWindFarm(GenericWindFarm):
    """
    A TrapWindFarm produces turbines on a grid but shape as a trapezoid. The params.yaml file determines
    how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                #                      # Description                | Units
                type: trap
                num_x:               # Number of columns             | -
                num_y:               # Number of rows                | -
                x_spacing: 300       # Space between columns         | m
                y_spacing: 300       # Space between row             | m
                expand_factor: 2     # How much the farm grows (>1) or shrinks (<1) between the first and last columns
                leading_center:[0,0] # location of the center of the leading column of turbines

                
        This will produce a 6x6 grid of turbines where the first column is 
        300 meters between turbines, and the last column is 600 m between 
        turbines. The first column will be centered on the point 0,0

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self,dom):

        self.name = "Trap Farm"
        super(TrapWindFarm, self).__init__(dom)

    def load_parameters(self):

        # Load from the yaml file
        self.num_x          = self.params["wind_farm"]["num_x"]
        self.num_y          = self.params["wind_farm"]["num_y"]
        self.x_spacing      = self.params["wind_farm"]["x_spacing"]
        self.y_spacing      = self.params["wind_farm"]["y_spacing"]
        self.expand_factor  = self.params["wind_farm"]["expand_factor"]
        self.leading_center = self.params["wind_farm"]["leading_center"]

    def compute_parameters(self):

        # compute number of turbines
        self.numturbs  = self.num_x*self.num_y

    def initialize_turbine_locations(self):

        ### Create the x and y coords ###
        if self.num_x > 1:
            x_list = np.linspace(0,1,self.num_x)
        else:
            raise ValueError("the number of column for a TrapWindFarm needs to be greater than 1")
        if self.num_y > 1:
            y_list = np.linspace(-1,1,self.num_y)
        else:
            raise ValueError("the number of rows for a TrapWindFarm needs to be greater than 1")

        ### Use the x and y coords to make a mesh grid ###
        X, Y = np.meshgrid(x_list,y_list)
        
        # Distort grid such that the first col remains unchanged but the last one has the right spacing
        Y = Y*(self.expand_factor*X-X+1)

        # Scale col to real locations
        X = self.x_spacing*(self.num_x-1)*(X)+self.leading_center[0]

        # Scale row to real locations
        Y = self.y_spacing*(self.num_y-1)*(Y/2)+self.leading_center[1]

        # Convert to lists
        x = X.flatten()
        y = Y.flatten()
            
        # Recompute the extents based on grid manipulations
        self.ex_x = [np.min(x), np.max(x)]
        self.ex_y = [np.min(y), np.max(y)]

        # Output some useful data
        self.fprint("Grid Size: {:d} x {:d}".format(self.num_y,self.num_x))
        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0],self.ex_x[1]))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0],self.ex_y[1]))

        return np.array([x,y]).T