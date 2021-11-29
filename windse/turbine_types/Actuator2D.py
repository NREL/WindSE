from . import GenericTurbine

class Actuator2D(GenericTurbine):

    def __init__(self,dom):
        super(ActuatorDisk, self).__init__(dom)

        ### special stuff here ###

    def LoadParameters():
        pass
        
    def UpdateControls():
        pass

    def Force():
        pass

    def ForceGradient():
        pass

    def Power():
        pass

    def PowerGradient():
        pass
