class Turbines(object):
    def __init__(self, x, y, HH, yaw, a, RD, W, ground):
        self.loc = loc # Location in Space
        self.RD  = RD  # Rotor Diameter
        self.W   = W   # Width of influence


from scipy.special import gamma
from mpmath import hyper

T_norm = 2.0*gamma(7.0/6.0)
D_norm = 2.0*gamma(7.0/6.0)
S_norm = (2.0+pi)/(2.0*pi)

h1 = float(hyper([],[1/3, 2/3, 5/6, 7/6],-(np.pi**6/46656)))
h2 = float(hyper([],[2/3, 7/6, 4/3, 3/2],-(np.pi**6/46656)))
h3 = float(hyper([],[4/3, 3/2, 5/3, 11/6],-(np.pi**6/46656)))
SD_norm = (1/3*np.pi**(3/2)*h1 - 1/15*np.pi**3*gamma(11/6)*h2 + gamma(7/6)*(1 + 1/360*np.pi**5*h3))/(.81831)
volNormalization = T_norm*SD_norm*W*R**(self.problem.dom.dim)