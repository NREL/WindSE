# import windse function
from . import GenericTurbine

# import dolfin functions
from . import Function, Constant

# import other modules
import numpy as np

class DisabledTurbine(GenericTurbine):

    def __init__(self, i,x,y,dom,imported_params=None):
        # Define the acceptable column of an wind_farm.csv imput file        
        self.yaml_inputs = ["HH", "RD", "yaw"]

        # Init turbine
        super(DisabledTurbine, self).__init__(i,x,y,dom,imported_params)
    

    def load_parameters(self):
        self.HH        = self.params["turbines"]["HH"]
        self.RD        = self.params["turbines"]["RD"]
        self.yaw       = self.params["turbines"]["yaw"]
    
    def create_controls(self):
        self.mx     = Constant(self.x, name="x_{:d}".format(self.index))
        self.my     = Constant(self.y, name="y_{:d}".format(self.index))
        self.myaw   = Constant(self.yaw, name="yaw_{:d}".format(self.index))
        # The control list is very important for extracting information from the optimzation step.
        self.controls_list = [] 

    def turbine_force(self,u,inflow_angle,fs):
        return Function(fs.V)

    def force_gradient(self):
        pass

    def power(self, u, inflow_angle):
        return 0.0

    def power_gradient(self):
        pass

    def prepare_saved_functions(self, func_list):
        return func_list