from . import GenericWindFarm
import numpy as np

class GridWindFarm(GenericWindFarm):
    """
    A GridWindFarm produces turbines on a grid. The params.yaml file determines
    how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                #                     # Description              | Units
                type: grid
                ex_x: [-1500, 1500]   # x-extent of the farm     | m
                ex_y: [-1500, 1500]   # y-extent of the farm     | m
                grid_rows: 6          # Number of rows           | -
                grid_cols: 6          # Number of columns        | -
                x_spacing: None       # x spacing between turbines in meters
                y_spacing: None       # y spacing between turbines in meters
                x_shear:   None       # x offset between rows in meters
                y_shear:   None       # y offset between columns in meters
                jitter:    0.0         # magnitude of random noise added to a gridded wind farm
        This will produce a 6x6 grid of turbines equally spaced within the 
        region [-1500, 1500]x[-1500, 1500].

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self,dom):

        self.name = "Grid Farm"
        super(GridWindFarm, self).__init__(dom)

    def load_parameters(self):

        # Load from the yaml file
        self.ex_x      = self.params["wind_farm"]["ex_x"]
        self.ex_y      = self.params["wind_farm"]["ex_y"]
        self.grid_rows = self.params["wind_farm"]["grid_rows"]
        self.grid_cols = self.params["wind_farm"]["grid_cols"]
        self.x_spacing = self.params["wind_farm"]["x_spacing"]
        self.y_spacing = self.params["wind_farm"]["y_spacing"]
        self.x_shear   = self.params["wind_farm"]["x_shear"]
        self.y_shear   = self.params["wind_farm"]["y_shear"]
        self.jitter    = self.params["wind_farm"]["jitter"]
        self.seed      = self.params["wind_farm"]["seed"]
        self.radius    = self.params["turbines"]["RD"]/2.0 # TODO: Should we get rid of the farm padding?
        self.numturbs  = self.grid_rows*self.grid_cols

        # Need to compute extents if only spacings are provided
        if self.ex_x is None:
            x_dist = self.x_spacing * (self.grid_cols - 1)
            y_dist = self.y_spacing * (self.grid_rows - 1)
            self.ex_x = [0., x_dist + 2 * self.radius]
            self.ex_y = [-y_dist / 2 - self.radius, y_dist / 2 + self.radius]
        else:
            self.ex_x   = [self.ex_x[0] + self.radius, self.ex_x[1] - self.radius]
            self.ex_y   = [self.ex_y[0] + self.radius, self.ex_y[1] - self.radius]

        # Output some useful data
        self.fprint("Grid Size: {:d} x {:d}".format(self.grid_rows,self.grid_cols))
        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0],self.ex_x[1]))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0],self.ex_y[1]))
        if self.jitter > 0.0:
            self.fprint("Amount of Jitter:   {: 1.2f}".format(self.jitter))
            self.fprint("Random Seed: " + repr(self.seed))

    def initialize_turbine_locations(self):
        ### Check if random seed is set ###
        if self.seed is not None:
            np.random.seed(self.seed)

        ### Create the x and y coords ###
        if self.grid_cols > 1:
            grid_x = np.linspace(self.ex_x[0],self.ex_x[1],self.grid_cols)
        else:
            grid_x = (self.ex_x[0]+self.ex_x[1])/2.0
        if self.grid_rows > 1:
            grid_y = np.linspace(self.ex_y[0],self.ex_y[1],self.grid_rows)
        else:
            grid_y = (self.ex_y[0]+self.ex_y[1])/2.0

        ### Use the x and y coords to make a mesh grid ###
        x, y = np.meshgrid(grid_x,grid_y)
        
        # Apply y shear if included in user yaml. Shear is in meters.
        if self.y_shear is not None:
            for idx in range(self.grid_cols):
                y[:, idx] += self.y_shear * idx
        
        # Apply x shear if included in user yaml. Shear is in meters.
        if self.x_shear is not None:
            for idx in range(self.grid_rows):
                x[idx, :] += self.x_shear * idx

        x = x.flatten()
        y = y.flatten()

        ### Apply Jitter ###
        if self.jitter > 0.0:
            if self.seed is not None:
                np.random.seed(self.seed)
            x += np.random.randn(self.numturbs)*self.jitter
            y += np.random.randn(self.numturbs)*self.jitter
            
        return np.array([x,y]).T