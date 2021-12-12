from . import GenericWindFarm

class ImportedWindFarm(GenericWindFarm):

    def __init__(self,dom):

        self.name = "Imported Farm"
        super(ImportedWindFarm, self).__init__(dom)

        ### special stuff here ###

    def load_parameters():
        pass

    def initialize_turbine_locations():
        pass