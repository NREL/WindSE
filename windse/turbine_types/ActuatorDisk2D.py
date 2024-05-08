# import windse function
from . import ActuatorDisk

# import dolfin functions
from . import SpatialCoordinate, as_vector, cos, sin, exp, sqrt, dot, assemble, dx
import numpy as np

class ActuatorDisk2D(ActuatorDisk):

    def __init__(self, i,x,y,dom,imported_params=None):
        super(ActuatorDisk2D, self).__init__(i,x,y,dom,imported_params)

    def power(self, u, inflow_angle):
        x=SpatialCoordinate(self.dom.mesh)

        mx = self.mx
        my = self.my
        mz = self.mz
        x0 = [mx,my,mz]
        W = self.thickness*1.0
        R = self.RD/2.0 
        ma = self.maxial
        yaw = self.myaw+self.inflow_angle
        A = np.pi*R**2.0
        C_tprime = 4*ma/(1-ma)
        C_pprime = 0.45/(1-ma)**3
        
        ### Rotate and Shift the Turbine ###
        xs = self.yaw_turbine(x,x0,yaw)
        u_d = u[0]*cos(yaw) + u[1]*sin(yaw)

        ### Create the function that represents the Thickness of the turbine ###
        T = exp(-pow((xs[0]/W),6.0))

        # WTGbase = Expression(("cos(yaw)","sin(yaw)"),yaw=yaw,degree=1)
        WTGbase = as_vector((cos(yaw),sin(yaw)))

        ### Create the function that represents the Disk of the turbine
        D = exp(-pow((pow((xs[1]/R),2)),6.0))

        volNormalization = assemble(T*D*dx)

        ### Create the function that represents the force ###
        if self.force == "constant":
            F = R*C_tprime    
        elif self.force == "sine":
            r = sqrt(xs[1]**2.0+xs[2]**2)
            F = R*C_tprime*(r/R*sin(np.pi*r/R)+0.5)/(.81831)

        J = (assemble(((0.5*A*C_pprime)**(1/3))*F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
        return J