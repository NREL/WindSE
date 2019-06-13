"""
The FunctionSpaceManager contains all the different types of function
spaces required for solve multiple classes of problems.
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import time

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

class GenericFunctionSpace(object):
    def __init__(self,dom):
        self.params = windse_parameters
        self.fprint = self.params.fprint
        self.dim = dom.dim

    def SetupSubspaces(self):
        self.V = self.W.sub(0).collapse()
        self.Q = self.W.sub(1).collapse()
        self.V0 = self.V.sub(0).collapse() 
        self.V1 = self.V.sub(1).collapse()
        if self.dim == 3:
            self.V2 = self.V.sub(2).collapse()

        if self.dim == 3: 
            self.VelocityAssigner = FunctionAssigner(self.V,[self.V0,self.V1,self.V2])
        else:
            self.VelocityAssigner = FunctionAssigner(self.V,[self.V0,self.V1])

        self.SolutionAssigner = FunctionAssigner(self.W,[self.V,self.Q])



class LinearFunctionSpace(GenericFunctionSpace):
    """
    The LinearFunctionSpace is made up of a vector function space for velocity
    and a scaler space for pressure. Both spaces are "CG1" or Linear Lagrange elements.
    """
    def __init__(self,dom):
        super(LinearFunctionSpace, self).__init__(dom)

        ### Create the function space ###
        fs_start = time.time()
        self.fprint("Creating Function Space",special="header")

        V = VectorElement('Lagrange', dom.mesh.ufl_cell(), 1) 
        Q = FiniteElement('Lagrange', dom.mesh.ufl_cell(), 1)
        self.W = FunctionSpace(dom.mesh, MixedElement([V,Q]))

        self.SetupSubspaces()

        self.fprint("Velocity DOFS: {:d}".format(self.V.dim()))
        self.fprint("Pressure DOFS: {:d}".format(self.Q.dim()))
        self.fprint("Total DOFS:    {:d}".format(self.W.dim()))
        
        fs_stop = time.time()
        self.fprint("Function Spaces Created: {:1.2f} s".format(fs_stop-fs_start),special="footer")

class TaylorHoodFunctionSpace(GenericFunctionSpace):
    """
    The TaylorHoodFunctionSpace is made up of a vector function space for velocity
    and a scalar space for pressure. The velocity function space is piecewise quadratic
    and the pressure function space is piecewise linear.
    """
    def __init__(self,dom):
        super(TaylorHoodFunctionSpace, self).__init__(dom)

        ### Create the function space ###
        fs_start = time.time()
        self.fprint("Creating Function Space",special="header")
        V = VectorElement('Lagrange', dom.mesh.ufl_cell(), 2) 
        Q = FiniteElement('Lagrange', dom.mesh.ufl_cell(), 1)
        self.W = FunctionSpace(dom.mesh, MixedElement([V,Q]))

        self.SetupSubspaces()

        self.fprint("Velocity DOFS: {:d}".format(self.V.dim()))
        self.fprint("Pressure DOFS: {:d}".format(self.Q.dim()))
        self.fprint("Total DOFS:    {:d}".format(self.W.dim()))

        fs_stop = time.time()
        self.fprint("Function Spaces Created: {:1.2f} s".format(fs_stop-fs_start),special="footer")

