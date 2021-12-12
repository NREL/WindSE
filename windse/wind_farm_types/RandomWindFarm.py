from . import GenericWindFarm
import numpy as np

class RandomWindFarm(GenericWindFarm):
    """
    A RandomWindFarm produces turbines located randomly with a defined 
    range. The params.yaml file determines how this grid is set up.

    Example:
        In the .yaml file you need to define::

            wind_farm: 
                #                     # Description              | Units
                type: random
                ex_x: [-1500, 1500]   # x-extent of the farm     | m
                ex_y: [-1500, 1500]   # y-extent of the farm     | m
                numturbs: 36          # Number of Turbines       | -
                seed: 15              # Random Seed for Numpy    | -

        This will produce a 36 turbines randomly located within the 
        region [-1500, 1500]x[-1500, 1500]. The seed is optional but 
        useful for reproducing test.

    Args: 
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """

    def __init__(self,dom):

        self.name = "Random Farm"
        super(RandomWindFarm, self).__init__(dom)

    def load_parameters(self):
        self.numturbs     = self.params["wind_farm"]["numturbs"]
        self.seed         = self.params["wind_farm"]["seed"]
        self.ex_x         = self.params["wind_farm"]["ex_x"]
        self.ex_y         = self.params["wind_farm"]["ex_y"]
        self.radius       = self.params["turbines"]["RD"]/2.0
        self.min_sep_dist = self.params["wind_farm"]["min_sep_dist"]*self.radius*2

    def initialize_turbine_locations(self):
        ### Check if random seed is set ###
        if self.seed is not None:
            np.random.seed(self.seed)

        ### Construct the x,y pairs using latin hypercubes
        rand_locations = self.build_random_samples(self.numturbs, self.ex_x, self.ex_y, self.min_sep_dist)

        self.fprint("X Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_x[0],self.ex_x[1]))
        self.fprint("Y Range: [{: 1.2f}, {: 1.2f}]".format(self.ex_y[0],self.ex_y[1]))
        self.fprint("Random Seed: " + repr(self.seed))

        return rand_locations

    def generate_random_point(self, x_range, y_range):
        
        rand_pt = np.zeros(2)
        
        # This assigns numbers in the range (x_range[0], x_range[1])
        rand_pt[0] = np.random.uniform(x_range[0], x_range[1])
        rand_pt[1] = np.random.uniform(y_range[0], y_range[1])

        return rand_pt

    def build_random_samples(self, N, x_range, y_range, min_dist, x_padding=None, y_padding=None):
        rand_samples = np.zeros((N, 2))

        if x_padding is not None:
            x_range[0] = x_range[0] + x_padding[0]
            x_range[1] = x_range[1] - x_padding[1]

        if y_padding is not None:
            y_range[0] = y_range[0] + y_padding[0]
            y_range[1] = y_range[1] - y_padding[1]
        
        # Declare the maximum number of attempts at placing a turbine          
        maximum_iterations = 50000

        for k in range(N):
            if k == 0:
                # The first turbine can always be added (guaranteed collision free)
                new_pt = self.generate_random_point(x_range, y_range)
                rand_samples[0, :] = new_pt

            else:
                # Additional turbines must be tested to enforce the minimum separation distance
                collision = True
                attempt = 0
                
                while collision == True:
                    new_pt = self.generate_random_point(x_range, y_range)
                    attempt += 1
                                    
                    dx_2 = (rand_samples[0:k, :] - new_pt)**2
                    dist_2 = np.sum(dx_2, axis = 1)
                                    
                    if np.amin(dist_2) < min_dist**2: # I think this need to be converted 
                        collision = True
                    else:
                        collision = False
                        rand_samples[k, :] = new_pt

                    if attempt > maximum_iterations:
                        # If the total numer of turbines couldn't be placed, raise an error or return the incomplete list
                        # (since the list is incomplete, numturbs needs to be updated with the reduced value)
                        # raise ValueError("Couldn't place point %d of %d after %d iterations." % (k+1, N, maximum_iterations))
                        self.fprint("WARNING: Couldn't place point %d of %d after %d iterations." % (k+1, N, maximum_iterations))
                        self.fprint("WARNING: Consider reducing the number of turbines or decreasing the minimum separation distance.")
                        self.fprint("WARNING: Proceeding with incomplete random farm, numturbs = %d turbines." % (k))
                        self.numturbs = k
                        return rand_samples[0:k, :]

        return rand_samples
