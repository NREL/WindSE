from . import GenericTurbine

class ActuatorLine(GenericTurbine):

    def __init__(self,dom):
        super(ActuatorLine, self).__init__(dom)

        ### special stuff here ###


    def get_baseline_chord(self):
        '''
        This function need to return the baseline chord as a numpy array with length num_blade_segments
        '''
        pass


    def load_parameters(self):

        self.max_chord = self.params["turbines"]["max_chord"]

        pass
        
    def create_controls(self):
        self.controls_list = ["x","y","yaw","chord","lift","drag"] # this is just part of the function as an example of the types of controls 
        pass

    def turbine_force(self,u,inflow_angle,fs,simTime):
        pass

    def power(self):
        pass

    def prepare_saved_functions(self, func_list):
        pass