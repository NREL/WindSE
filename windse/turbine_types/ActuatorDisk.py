from . import GenericTurbine
from . import Constant

class ActuatorDisk(GenericTurbine):

    def __init__(self, i,x,y,dom):
        super(ActuatorDisk, self).__init__(i,x,y,dom)

        ### special stuff here ###

    def LoadParameters(self):
        self.HH        = self.params["turbines"]["HH"]
        self.RD        = self.params["turbines"]["RD"]
        self.thickness = self.params["turbines"]["thickness"]
        self.yaw       = self.params["turbines"]["yaw"]
        self.axial     = self.params["turbines"]["axial"]
    
    def CreateControls(self):
        self.mx   = Constant(self.x, name="x_{:d}".format(self.index))
        self.my   = Constant(self.y, name="y_{:d}".format(self.index))
        self.myaw = Constant(self.yaw, name="yaw_{:d}".format(self.index))
        self.ma   = Constant(self.axial, name="axial_{:d}".format(self.index))

    def UpdateControls(self):
        self.mx.assign(self.x) 
        self.my.assign(self.y)    
        self.myaw.assign(self.yaw)  
        self.ma.assign(self.axial)    

    def Force(self):
        pass

    def ForceGradient(self):
        pass

    def Power(self):
        pass

    def PowerGradient(self):
        pass
