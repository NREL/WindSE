from . import GenericWindFarm

class EmptyWindFarm(GenericWindFarm):

    def __init__(self,dom):

        self.name = "Empty Farm"
        super(EmptyWindFarm, self).__init__(dom)

    def LoadParameters():
        pass

    def InitializeTurbineLocations():
        pass