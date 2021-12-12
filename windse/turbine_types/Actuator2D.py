from . import GenericTurbine

class Actuator2D(GenericTurbine):

    def __init__(self,dom):
        super(ActuatorDisk, self).__init__(dom)

        ### special stuff here ###

    def load_parameters(self):
        pass
        
    def create_controls(self):
        pass

    def update_controls(self):
        pass  

    def turbine_force(self,u,inflow_angle,fs,simTime):
        pass

    def power(self):
        pass
