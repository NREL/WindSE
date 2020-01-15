"""
The ProblemManager contains all of the 
different classes of problems that windse can solve
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

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
    def __init__(self,domain,windfarm,function_space,boundary_data):
        ### save a reference of option and create local version specifically of domain options ###
        self.params = windse_parameters
        self.dom  = domain
        self.farm = windfarm
        self.fs   = function_space 
        self.bd  = boundary_data
        self.tf_first_save = True
        self.fprint = self.params.fprint

        ### Append the farm ###
        self.params.full_farm = self.farm

        ### add some helper items for dolfin_adjoint_helper.py ###
        self.params.ground_fx = self.dom.Ground
        self.params.full_hh = self.farm.HH

    def ComputeTurbineForce(self,theta):

        ### Compute the relative yaw angle ###
        if theta is not None:
            inflow_angle = theta-self.dom.init_wind
        else:
            inflow_angle = 0.0

        ### Create the turbine force function ###
        if self.farm.turbine_method == "dolfin":
                self.tf1, self.tf2, self.tf3 = self.farm.DolfinTurbineForce(self.fs,self.dom.mesh,inflow_angle=inflow_angle)

        elif self.farm.turbine_method == "numpy":
                self.tf1, self.tf2, self.tf3 = self.farm.NumpyTurbineForce(self.fs,self.dom.mesh,inflow_angle=inflow_angle)
        else:
            raise ValueError("Unknown turbine method: "+self.farm.turbine_method)
        
        ### Convolve TF with u ###
        tf = self.tf1*self.u_next[0]**2+self.tf2*self.u_next[1]**2+self.tf3*self.u_next[0]*self.u_next[1]

        return tf

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

        adj_stop = time.time()
        self.fprint("Wind Angle Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")

        # self.tf.assign(tf_temp)

    def ChangeWindSpeed(self,speed):
        adj_start = time.time()
        self.fprint("Adjusting Wind Speed",special="header")
        self.fprint("New Speed: {:1.8f} m/s".format(speed/self.dom.xscale))
        self.bd.HH_vel = speed
        adj_stop = time.time()
        self.fprint("Wind Speed Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")

class StabilizedProblem(GenericProblem):
    """
    The StabilizedProblem setup everything required for solving Navier-Stokes with 
    a stabilization term

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_conditions):
        super(StabilizedProblem, self).__init__(domain,windfarm,function_space,boundary_conditions)
        
        ### Create Functional ###
        self.ComputeFunctional()


    def ComputeFunctional(self,theta=None):
        self.fprint("Setting Up Stabilized Problem",special="header")

        ### Create the test/trial/functions ###
        self.up_next = Function(self.fs.W)
        self.u_next,self.p_next = split(self.up_next)
        v,q = TestFunctions(self.fs.W)

        ### Set the x scaling ###
        Sx = self.dom.xscale

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_next.assign(self.bd.u0)

        ### Create the turbine force ###
        self.tf = self.ComputeTurbineForce(theta)

        ### These constants will be moved into the params file ###
        nu = self.params["problem"].get("viscosity",.1)
        f = Constant((0.0,)*self.dom.dim)
        f.rename("f","f")
        
        vonKarman=0.41
        lmax=15
        mlDenom = 8.
        eps=Constant(1.0)
        eps.rename("eps","eps")

        self.fprint("Viscosity:                 {:1.2e}".format(float(nu)))
        self.fprint("Mixing Length Scale:       {:1.2e}".format(float(mlDenom)))
        self.fprint("Stabilization Coefficient: {:1.2e}".format(float(eps)))

        ### Calculate the stresses and viscosities ###
        S = sqrt(2*inner(0.5*(grad(self.u_next)+grad(self.u_next).T),0.5*(grad(self.u_next)+grad(self.u_next).T)))

        ### Create l_mix based on distance to the ground ###
        if self.dom.dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,(1.+np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,lmax)))
        else:
            l_mix = Constant((self.farm.HH[0]/Sx)/mlDenom)
        # l_mix = Expression("x[2]/8.",degree=2)
        l_mix.rename("l_mix","l_mix")
        
        ### Calculate nu_T
        self.nu_T=l_mix**2.*S

        print(assemble(dot(self.tf,self.u_next)*dx))
        # exit()


        ### Create the functional ###
        # if self.farm.yaw[0]**2 > 1e-4:
        #     self.F = inner(grad(self.u_next)*self.u_next, v)*dx + (nu+self.nu_T)*inner(grad(self.u_next), grad(v))*dx - inner(div(v),self.p_next)*dx - inner(div(self.u_next),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
        # else :
        # self.F = inner(grad(self.u_next)*self.u_next, v)*dx + Sx*Sx*inner(grad(self.u_next), grad(v))*dx - inner(div(v),self.p_next)*dx - inner(div(self.u_next),q)*dx - inner(f,v)*dx# + inner(self.tf,v)*dx 
        self.F = inner(grad(self.u_next)*self.u_next, v)*dx + Sx*Sx*(nu+self.nu_T)*inner(grad(self.u_next), grad(v))*dx - inner(div(v),self.p_next)*dx - inner(div(self.u_next),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
        # self.F_sans_tf =  (1.0)*inner(grad(self.u_next), grad(v))*dx - inner(div(v),self.p_next)*dx - inner(div(self.u_next),q)*dx - inner(f,v)*dx
        # self.F = inner(grad(self.u_next)*self.u_next, v)*dx + (nu+self.nu_T)*inner(grad(self.u_next), grad(v))*dx - inner(div(v),self.p_next)*dx - inner(div(self.u_next),q)*dx - inner(f,v)*dx + inner(self.tf*(self.u_next[0]**2+self.u_next[1]**2),v)*dx 

        ### Add in the Stabilizing term ###
        # stab = - eps*inner(grad(q), grad(self.p_next))*dx - eps*inner(grad(q), dot(grad(self.u_next), self.u_next))*dx 
        stab = - eps*inner(grad(q), grad(self.p_next))*dx - eps*inner(grad(q), dot(grad(self.u_next), self.u_next))*dx 
        # stab_sans_tf = - eps*inner(grad(q), grad(self.p_next))*dx 
        self.F += stab
        # self.F_sans_tf += stab
        self.fprint("Stabilized Problem Setup",special="footer")


class TaylorHoodProblem(GenericProblem):
    """
    The TaylorHoodProblem sets up everything required for solving Navier-Stokes 

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self,domain,windfarm,function_space,boundary_conditions):
        super(TaylorHoodProblem, self).__init__(domain,windfarm,function_space,boundary_conditions)

        ### Create Functional ###
        self.ComputeFunctional()

    def ComputeFunctional(self,theta=None):
        self.fprint("Setting Up Taylor-Hood Problem",special="header")

        ### These constants will be moved into the params file ###
        nu = self.params["problem"].get("viscosity",0.1)
        f = Constant((0.0,)*self.dom.dim)
        vonKarman=0.41
        lmax=15
        mlDenom = 8.
        eps=Constant(0.01)

        self.fprint("Viscosity:           {:1.2e}".format(float(nu)))
        self.fprint("Mixing Length Scale: {:1.2e}".format(float(mlDenom)))

        ### Create the test/trial/functions ###
        self.up_next = Function(self.fs.W)
        self.u_next,self.p_next = split(self.up_next)
        v,q = TestFunctions(self.fs.W)

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_next.assign(self.bd.u0)

        ### Calculate the stresses and viscosities ###
        S = sqrt(2.*inner(0.5*(grad(self.u_next)+grad(self.u_next).T),0.5*(grad(self.u_next)+grad(self.u_next).T)))

        ### Create l_mix based on distance to the ground ###
        if self.dom.dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:],(1.+np.divide(vonKarman*self.bd.depth.vector()[:],lmax)))
        else:
            l_mix = Constant(self.farm.HH[0]/mlDenom)

        ### Calculate nu_T
        self.nu_T=l_mix**2.*S

        ### Create the turbine force ###
        self.tf = self.ComputeTurbineForce(theta)

        ### Create the functional ###
        self.F = inner(grad(self.u_next)*self.u_next, v)*dx + (nu+self.nu_T)*inner(grad(self.u_next), grad(v))*dx - inner(div(v),self.p_next)*dx - inner(div(self.u_next),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
    
        self.fprint("Taylor-Hood Problem Setup",special="footer")

class UnsteadyProblem(GenericProblem):
    """
    The UnsteadyProblem sets up everything required for solving Navier-Stokes using
    a fractional-step method with an adaptive timestep size

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self, domain, windfarm, function_space, boundary_conditions):
        super(UnsteadyProblem, self).__init__(domain, windfarm, function_space, boundary_conditions)
        self.fprint("Setting Up Unsteady Problem", special="header")

        # ================================================================

        # Define fluid properties
        # FIXME: These should probably be set in params.yaml input filt
        mu = 1/10000
        rho = 1
        mu_c = Constant(mu)
        rho_c = Constant(rho)

        # Define time step size (this value is used only for step 1 if adaptive timestepping is used)
        # FIXME: change variable name to avoid confusion within dolfin adjoint
        self.dt = 0.1*self.dom.mesh.hmin()/self.bd.HH_vel
        self.dt_c  = Constant(self.dt)

        self.fprint("Viscosity: {:1.2e}".format(float(mu)))
        self.fprint("Density:   {:1.2e}".format(float(rho)))

        # Define trial and test functions for velocity
        u = TrialFunction(self.fs.V)
        v = TestFunction(self.fs.V)

        # Define trial and test functions for pressure
        p = TrialFunction(self.fs.Q)
        q = TestFunction(self.fs.Q)

        # Define functions for velocity solutions
        # >> _k = current (step k)
        # >> _k1 = previous (step k-1)
        # >> _k2 = double previous (step k-2)
        self.u_k = Function(self.fs.V)
        self.u_k1 = Function(self.fs.V)
        self.u_k2 = Function(self.fs.V)

        # Seed previous velocity fields with the chosen initial condition
        self.u_k1.assign(self.bd.bc_velocity)
        self.u_k2.assign(self.bd.bc_velocity)

        # Define functions for pressure solutions
        # >> _k = current (step k)
        # >> _k1 = previous (step k-1)
        self.p_k  = Function(self.fs.Q)
        self.p_k1 = Function(self.fs.Q)

        # Seed previous pressure fields with the chosen initial condition
        self.p_k1.assign(self.bd.bc_pressure)

        # ================================================================

        # Crank-Nicolson velocity
        U_CN  = 0.5*(u + self.u_k1)

        # Adams-Bashforth projected velocity
        U_AB = 1.5*self.u_k1 - 0.5*self.u_k2

        # ================================================================

        # Calculate eddy viscosity, if using
        use_eddy_viscosity = True

        if use_eddy_viscosity:
            # Define filter scale
            filter_scale = CellVolume(self.dom.mesh)**(1.0/self.dom.dim)

            # Strain rate tensor, 0.5*(du_i/dx_j + du_j/dx_i)
            Sij = sym(nabla_grad(U_AB))

            # sqrt(Sij*Sij)
            strainMag = (2.0*inner(Sij, Sij))**0.5

            # Smagorinsky constant, typically around 0.17
            Cs = 0.17

            # Eddy viscosity
            self.nu_T = Cs**2 * filter_scale**2 * strainMag
        else:
            self.nu_T = Constant(0)

        # ================================================================

        # FIXME: This up_next function is only present to avoid errors  
        # during assignments in GenericSolver.__init__

        # Create the combined function space
        self.up_next = Function(self.fs.W)

        # Create the turbine force
        # FIXME: Should this be set by a numpy array operation or a fenics function?
        # self.tf = self.farm.TurbineForce(self.fs, self.dom.mesh, self.u_k2)
        self.tf = Function(self.fs.V)

        # self.u_k2.vector()[:] = 0.0
        # self.u_k1.vector()[:] = 0.0


        # ================================================================

        # Define variational problem for step 1: tentative velocity
        F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
           + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
           + (mu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
           + dot(nabla_grad(self.p_k1), v)*dx \
           - dot(-self.tf, v)*dx

        self.a1 = lhs(F1)
        self.L1 = rhs(F1)

        # Define variational problem for step 2: pressure correction
        self.a2 = dot(nabla_grad(p), nabla_grad(q))*dx
        self.L2 = dot(nabla_grad(self.p_k1), nabla_grad(q))*dx - (1.0/self.dt_c)*div(self.u_k)*q*dx

        # Define variational problem for step 3: velocity update
        self.a3 = dot(u, v)*dx
        self.L3 = dot(self.u_k, v)*dx - self.dt_c*dot(nabla_grad(self.p_k - self.p_k1), v)*dx
    
        # ================================================================

        self.fprint("Unsteady Problem Setup",special="footer")
