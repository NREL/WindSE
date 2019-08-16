"""
The PhysicsManager contains different turbulence models.
"""

import __main__
import os

### Get the name of program importing this package ###
try:
 main_file = os.path.basename(__main__.__file__)
except:
 main_file = ""

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np
    import time

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *


class GenericTurbulenceModel(object):
    """
    """

    def __init__(self,):
        self.model = ''
        self.params = {}
        self.nu_T = None

    def calculate_nut(self, u_next, p_next, nu, fs, **kwargs):
        return self.nu_T


class MixingLengthTurbulenceModel(GenericTurbulenceModel):
    """
    """

    def __init__(self,):
        super(MixingLengthTurbulenceModel, self).__init__()
        self.model = 'mixing_length'
        self.params['vonKarman'] = 0.41
        self.params['lmax'] = 15
        self.params['mlDenom'] = 8.

    def calculate_nut(self, u_next, p_next, nu, fs, dim, HH, depth):
        windse_parameters.fprint("Mixing Length Scale:       {:1.2e}".format(float(self.params['mlDenom'])))

        ### Calculate the stresses and viscosities ###
        S = sqrt(2.*inner(0.5*(grad(u_next)+grad(u_next).T),0.5*(grad(u_next)+grad(u_next).T)))

        ### Create l_mix based on distance to the ground ###
        if dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(fs.Q)
            l_mix.vector()[:] = np.divide(self.params['vonKarman']*depth,(1.+np.divide(self.params['vonKarman']*depth,self.params['lmax'])))
        else:
            l_mix = Constant(HH/self.params['mlDenom'])

        ### Calculate nu_T
        return l_mix**2.*S


class SpecifiedNuT(GenericTurbulenceModel):
    """
    """

    def __init__(self,):
        super(SpecifiedNuT, self).__init__()
        self.model = 'specified_nuT'
        self.params['nuT_file'] = windse_parameters["physics"].get("nuT_file", "nu_T.xml")

    def calculate_nut(self, u_next, p_next, nu, fs):
        windse_parameters.fprint("Reading nuT file: {}".format(self.params['nuT_file']))
        return Function(fs.Q, self.params['nuT_file'])
