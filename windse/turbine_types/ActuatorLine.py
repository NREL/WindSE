from . import GenericTurbine

class ActuatorLine(GenericTurbine):

    def __init__(self,dom):
        super(ActuatorLine, self).__init__(dom)

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
