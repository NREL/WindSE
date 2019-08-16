"""
The ProblemManager contains all of the
different classes of problems that windse can solve
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


class GenericProblem(object):
    """
    A GenericProblem contains on the basic functions required by all problem objects.

    Args:
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_data,turbulence_model):
        ### Save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.dom  = domain
        self.farm = windfarm
        self.fs   = function_space
        self.bd  = boundary_data
        self.tm = turbulence_model
        # self.tf_first_save = True
        self.fprint = self.params.fprint

        ### Append the farm ###
        self.params.full_farm = self.farm

        ### Add some helper items for dolfin_adjoint_helper.py ###
        self.params.ground_fx = self.dom.Ground
        self.params.full_hh = self.farm.HH

        ### Fluid property ###
        self.nu = Constant(0.1)

        ### Turbulence model additional arguments ###
        self.tm_kwargs = {}
        if self.tm.model=='mixing_length':
            self.tm_kwargs['dim'] = self.dom.dim
            self.tm_kwargs['HH'] = self.farm.HH[0]
            if self.tm_kwargs['dim']==3:
                self.tm_kwargs['depth'] = self.bd.depth.vector()[:]
            else:
                self.tm_kwargs['depth'] = None
        elif self.tm.model=='specified_nuT':
            pass
        else:
            raise ValueError('Turbulence model "{}" not implemented.'.format(self.tm.model))

    def ChangeWindAngle(self,theta):
        """
        This function recomputes all necessary components for a new wind direction

        Args:
            theta (float): The new wind angle in radians
        """
        adj_start = time.time()
        self.fprint("Adjusting Wind Angle",special="header")
        self.fprint("New Angle: {:1.8f} rads".format(theta))
        self.dom.RecomputeBoundaryMarkers(theta)
        self.bd.RecomputeVelocity(theta)
        self.ComputeFunctional(theta)

        # self.tf.assign(tf_temp)
        adj_stop = time.time()
        self.fprint("Wind Angle Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")


class TaylorHoodProblem(GenericProblem):
    """
    The TaylorHoodProblem sets up everything required for solving Navier-Stokes

    Args:
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_conditions,turbulence_model):
        super(TaylorHoodProblem, self).__init__(domain,windfarm,function_space,boundary_conditions,turbulence_model)
        self.fprint("Setting Up Taylor-Hood Problem",special="header")

        ### Create Functional ###
        self.ComputeFunctional()

        self.fprint("Taylor-Hood Problem Setup",special="footer")

    def ComputeFunctional(self, theta=None):
        ### Create the test/trial/functions ###
        self.up_next = Function(self.fs.W)
        self.u_next,self.p_next = split(self.up_next)
        self.v,self.q = TestFunctions(self.fs.W)

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_next.assign(self.bd.u0)

        ### Calculate nu_T
        self.nu_T = self.tm.calculate_nut(self.u_next,self.p_next,self.nu, self.fs, **self.tm_kwargs)

        ### Create the turbine force ###
        if theta is not None:
            self.tf = self.farm.TurbineForce(self.fs,self.dom.mesh,self.u_next,delta_yaw=(theta-self.dom.init_wind))
        else:
            self.tf = self.farm.TurbineForce(self.fs,self.dom.mesh,self.u_next)

        ### Create the functional ###
        f = Constant((0.0,)*self.dom.dim)
        if self.farm.yaw[0]**2 > 1e-4:
            self.F = inner(grad(self.u_next)*self.u_next, self.v)*dx + (self.nu+self.nu_T)*inner(grad(self.u_next), grad(self.v))*dx - inner(div(self.v),self.p_next)*dx - inner(div(self.u_next),self.q)*dx - inner(f,self.v)*dx + inner(self.tf,self.v)*dx
        else :
            self.F = inner(grad(self.u_next)*self.u_next, self.v)*dx + (self.nu+self.nu_T)*inner(grad(self.u_next), grad(self.v))*dx - inner(div(self.v),self.p_next)*dx - inner(div(self.u_next),self.q)*dx - inner(f,self.v)*dx + inner(self.tf*(self.u_next[0]**2+self.u_next[1]**2),self.v)*dx


class StabilizedProblem(TaylorHoodProblem, GenericProblem):
    """
    The StabilizedProblem setup everything required for solving Navier-Stokes with
    a stabilization term

    Args:
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_conditions,turbulence_model):
        windse_parameters.fprint("Setting Up Stabilized Problem",special="header")

        super(StabilizedProblem, self).__init__(domain,windfarm,function_space,boundary_conditions,turbulence_model)

        self.fprint("Adding stabilization term")

        self.fprint("Stabilized Problem Setup",special="footer")


    def ComputeFunctional(self,theta=None):
        super().ComputeFunctional(theta)

        ### Add in the Stabilizing term ###
        eps = Constant(0.01)
        stab = - eps*inner(grad(self.q), grad(self.p_next))*dx - eps*inner(grad(self.q), dot(grad(self.u_next), self.u_next))*dx
        self.F += stab
