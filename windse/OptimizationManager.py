"""
The OptimizationManager submodule contains all the required function for
optimizing via dolfin-adjoint. To use dolfin-adjoin set::

    general:
        dolfin_adjoint: True

in the param.yaml file.

Todo:
    * Read through an update the docstrings for these functions.
    * Create specific optimization classes.
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
    import copy

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

class Optimizer(object):
    """
    A GenericProblem contains on the basic functions required by all problem objects.

    Args:
        dom (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
    """
    def __init__(self, solver):
        ### save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.solver = solver
        self.problem = solver.problem
        self.farm = solver.problem.farm
        self.fprint = self.params.fprint

        self.fprint("Setting Up Optimizer",special="header")
        self.fprint("Controls: {0}".format(self.params["optimization"]["controls"]))
        self.CreateControls()
        self.fprint("Number of Controls: {:d}".format(len(self.controls)))
        self.fprint("Define Bounds")
        self.CreateBounds()
        self.fprint("Define Optimizing Functional")
        self.PowerFunctional()
        self.fprint("Optimizer Setup",special="footer")

    def CreateControls(self):
        self.controls = []
        self.names = []
        self.indexes = [[],[],[],[]]
        self.init_vals = []
        j = 0
        if "layout" in self.params["optimization"]["controls"]:
            for i in range(self.farm.numturbs):
                self.indexes[0].append(j)
                j+=1
                self.names.append("x_"+repr(i))
                self.controls.append(Control(self.farm.mx[i]))
                self.init_vals.append(self.farm.mx[i])

                self.indexes[1].append(j)
                j+=1
                self.names.append("y_"+repr(i))
                self.controls.append(Control(self.farm.my[i]))
                self.init_vals.append(self.farm.my[i])

        if "yaw" in self.params["optimization"]["controls"]:
            for i in range(self.farm.numturbs):
                self.indexes[2].append(j)
                j+=1
                self.names.append("yaw_"+repr(i))
                self.controls.append(Control(self.farm.myaw[i]))
                self.init_vals.append(self.farm.myaw[i])

        if "axial" in self.params["optimization"]["controls"]:
            for i in range(self.farm.numturbs):
                self.indexes[3].append(j)
                j+=1
                self.names.append("axial_"+repr(i))
                self.controls.append(Control(self.farm.ma[i]))
                self.init_vals.append(self.farm.ma[i])

    def CreateBounds(self):
        lower_bounds = []
        upper_bounds = []

        if "layout" in self.params["optimization"]["controls"]:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant((self.farm.ex_x[0] + self.farm.radius[i])))
                lower_bounds.append(Constant((self.farm.ex_y[0] + self.farm.radius[i])))
                upper_bounds.append(Constant((self.farm.ex_x[1] - self.farm.radius[i])))
                upper_bounds.append(Constant((self.farm.ex_y[1] - self.farm.radius[i])))

        if "yaw" in self.params["optimization"]["controls"]:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(-pi/4.))
                upper_bounds.append(Constant(pi/4.))

        if "axial" in self.params["optimization"]["controls"]:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(1.))

        self.bounds = [lower_bounds,upper_bounds]

    def AssignControls(self):
        if "layout" in self.params["optimization"]["controls"]:
            self.farm.mx = self.controls[self.indexes[0]]
            self.farm.my = self.controls[self.indexes[1]]
        if "yaw" in self.params["optimization"]["controls"]:
            self.farm.myaw = self.controls[self.indexes[2]]
        if "axial" in self.params["optimization"]["controls"]:
            self.farm.ma = self.controls[self.indexes[3]]



    def PowerFunctional(self):
        """
        Creates the power functional that will be optimized

        Args:
            tf (dolfin.Function): Turbine Force function
            u (dolfin.Function): Velocity vector.
        """
        #how to handle rotation?
        # J=Functional(tf*u[0]**3*dx)
        self.J=assemble(-dot(self.problem.tf,self.solver.u_next)*dx)
        self.Jhat = ReducedFunctional(self.J, self.controls)

    def ListControls(self,m):
        for i,val in enumerate(m):
            self.fprint(self.names[i]+": "+repr(m[i]))

    def Optimize(self):

        self.fprint("Beginning Optimization",special="header")
        m_opt=minimize(self.Jhat, method="L-BFGS-B", options = {"disp": True}, bounds = self.bounds, callback = self.ListControls)

        self.fprint("Assigning New Values")
        self.AssignControls()

        self.fprint("Solving With New Values")
        self.solver.Solve()

        self.fprint("Optimization Finished",special="header")

        return m_opt

    def TaylorTest(self):

        self.fprint("Beginning Taylor Test",special="header")

        h = [Constant(0.001)]*(len(self.controls))

        conv_rate = taylor_test(self.Jhat, self.init_vals, h)

        self.fprint("Convergence Rates:")
        self.fprint("")
        self.fprint(conv_rate)
        self.fprint("")

        self.fprint("Taylor Test Finished",special="footer")

        return conv_rate
