from . import GenericWindFarm

class ImportedWindFarm(GenericWindFarm):

    def __init__(self,dom):

        self.name = "Imported Farm"
        super(ImportedWindFarm, self).__init__(dom)

        ### special stuff here ###

    def LoadParameters():
        pass

    def InitializeTurbineLocations():
        pass