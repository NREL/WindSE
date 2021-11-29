from . import GenericWindFarm


class GridWindFarm(GenericWindFarm):

    def __init__(self,dom):

        self.name = "Grid Farm"
        super(GridWindFarm, self).__init__(dom)

        ### special stuff here ###

    def LoadParameters():
        pass

    def InitializeTurbineLocations():
        pass