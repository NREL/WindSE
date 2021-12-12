from . import GenericWindFarm

class EmptyWindFarm(GenericWindFarm):

    def __init__(self,dom):

        self.name = "Empty Farm"
        super(EmptyWindFarm, self).__init__(dom)

    def load_parameters():
        pass

    def initialize_turbine_locations():
        pass