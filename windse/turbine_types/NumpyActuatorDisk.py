from . import ActuatorDisk
import numpy as np

class NumpyActuatorDisk(ActuatorDisk):

    def __init__(self, i,x,y,dom):
        super(NumpyActuatorDisk, self).__init__(i,x,y,dom)

    def build_actuator_disk(self,inflow_angle):
        raise NotImplementedError(type(self))
        

    def turbine_force(self,u,inflow_angle,simTime=0):
        raise NotImplementedError(type(self))
