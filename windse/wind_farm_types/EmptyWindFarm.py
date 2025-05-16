from . import GenericWindFarm

class EmptyWindFarm(GenericWindFarm):

    def __init__(self):

        self.name = "Empty Farm"
        super(EmptyWindFarm, self).__init__()

    def load_parameters(self):
        pass

    def compute_parameters(self):
        self.numturbs = 0
        pass

    def initialize_turbine_locations(self):
        return []