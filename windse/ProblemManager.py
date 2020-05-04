"""
The ProblemManager contains all of the 
different classes of problems that windse can solve
"""

import __main__
import os

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np
    import time

    ### Import the cumulative parameters ###
    from windse import windse_parameters, CalculateActuatorLineTurbineForces
    # from memory_profiler import memory_usage

    ### Check if we need dolfin_adjoint ###
    if windse_parameters.dolfin_adjoint:
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

        ### Update attributes based on params file ###
        for key, value in self.params["problem"].items():
            setattr(self,key,value)

    def ComputeTurbineForce(self,u,inflow_angle,simTime=0.0):

        ### Compute the relative yaw angle ###
        if inflow_angle is not None:
            inflow_angle = inflow_angle-self.dom.inflow_angle
        else:
            inflow_angle = 0.0

        ### Create the turbine force function ###
        if self.farm.turbine_method == "dolfin":
            self.tf1, self.tf2, self.tf3 = self.farm.DolfinTurbineForce(self.fs,self.dom.mesh,inflow_angle=inflow_angle)

        elif self.farm.turbine_method == "numpy":
            self.tf1, self.tf2, self.tf3 = self.farm.NumpyTurbineForce(self.fs,self.dom.mesh,inflow_angle=inflow_angle)

        elif self.farm.turbine_method == 'alm':
            self.rpm = self.params["wind_farm"]["rpm"]
            self.num_blade_segments = 10
            self.mcl = []
            self.mcd = []

            # cl = np.linspace(0.0, 2.0, self.num_blade_segments)
            # cd = np.linspace(2.0, 0.0, self.num_blade_segments)
            cl = np.ones(self.num_blade_segments)
            cd = np.ones(self.num_blade_segments)

            for k in range(self.num_blade_segments):
                self.mcl.append(Constant(cl[k]))
                self.mcd.append(Constant(cd[k]))
            tf = CalculateActuatorLineTurbineForces(self, simTime)
        else:
            raise ValueError("Unknown turbine method: "+self.farm.turbine_method)
        
        ### Convolve TF with u ###
        if self.farm.turbine_method != 'alm':
            tf = self.tf1*u[0]**2+self.tf2*u[1]**2+self.tf3*u[0]*u[1]

        return tf

    def ChangeWindAngle(self,inflow_angle):
        """
        This function recomputes all necessary components for a new wind direction

        Args: 
            inflow_angle (float): The new wind angle in radians
        """
        adj_start = time.time()
        self.fprint("Adjusting Wind Angle",special="header")
        self.fprint("New Angle: {:1.8f} rads".format(inflow_angle))
        self.dom.RecomputeBoundaryMarkers(inflow_angle)
        self.bd.RecomputeVelocity(inflow_angle)
        self.ComputeFunctional(inflow_angle)

        adj_stop = time.time()
        self.fprint("Wind Angle Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")

        # self.tf.assign(tf_temp)

    def ChangeWindSpeed(self,inflow_speed):
        adj_start = time.time()
        self.fprint("Adjusting Wind Speed",special="header")
        self.fprint("New Speed: {:1.8f} m/s".format(inflow_speed/self.dom.xscale))
        self.bd.HH_vel = inflow_speed
        adj_stop = time.time()
        self.fprint("Wind Speed Adjusted: {:1.2f} s".format(adj_stop-adj_start),special="footer")

    def UpdateActuatorLineControls(self, c_lift = None, c_drag = None):

        if c_lift is not None:
            cl = np.array(c_lift, dtype = float)
        if c_drag is not None:
            cd = np.array(c_drag, dtype = float)

        for k in range(self.num_blade_segments):
            self.mcl[k] = Constant(cl[k])
            self.mcd[k] = Constant(cd[k])


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


    def ComputeFunctional(self,inflow_angle=None):
        self.fprint("Setting Up Stabilized Problem",special="header")

        ### Create the test/trial/functions ###
        self.up_k = Function(self.fs.W)
        self.u_k,self.p_k = split(self.up_k)
        v,q = TestFunctions(self.fs.W)

        ### Set the x scaling ###
        Sx = self.dom.xscale

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_k.assign(self.bd.u0)

        # mem0=memory_usage()[0]
        # mem_out, self.tf = memory_usage((self.ComputeTurbineForce,(self.u_k,inflow_angle),{}),max_usage=True,retval=True,max_iterations=1)
        # self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))
        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)

        ### These constants will be moved into the params file ###
        f = Constant((0.0,)*self.dom.dim)
        f.rename("f","f")
        
        nu = self.viscosity
        vonKarman=0.41
        eps=Constant(1.0)
        eps.rename("eps","eps")

        self.fprint("Viscosity:                 {:1.2e}".format(float(self.viscosity)))
        self.fprint("Max Mixing Length:         {:1.2e}".format(float(self.lmax)))
        self.fprint("Stabilization Coefficient: {:1.2e}".format(float(eps)))

        ### Calculate the stresses and viscosities ###
        S = sqrt(2*inner(0.5*(grad(self.u_k)+grad(self.u_k).T),0.5*(grad(self.u_k)+grad(self.u_k).T)))

        ### Create l_mix based on distance to the ground ###
        if self.dom.dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,(1.+np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,self.lmax)))
        else:
            l_mix = Constant(vonKarman*self.farm.HH[0]/(1+(vonKarman*self.farm.HH[0]/self.lmax)))
        # l_mix = Expression("x[2]/8.",degree=2)
        l_mix.rename("l_mix","l_mix")
        
        ### Calculate nu_T
        self.nu_T=l_mix**2.*S

        ### Create the functional ###
        # if self.farm.yaw[0]**2 > 1e-4:
        #     self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
        # else :
        # self.F = inner(grad(self.u_k)*self.u_k, v)*dx + Sx*Sx*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx# + inner(self.tf,v)*dx 
        self.F = inner(grad(self.u_k)*self.u_k, v)*dx + Sx*Sx*(nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 
        # self.F_sans_tf =  (1.0)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx
        # self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf*(self.u_k[0]**2+self.u_k[1]**2),v)*dx 

        ### Add in the Stabilizing term ###
        # stab = - eps*inner(grad(q), grad(self.p_k))*dx - eps*inner(grad(q), dot(grad(self.u_k), self.u_k))*dx 
        stab = - eps*inner(grad(q), grad(self.p_k))*dx - eps*inner(grad(q), dot(grad(self.u_k), self.u_k))*dx 
        # stab_sans_tf = - eps*inner(grad(q), grad(self.p_k))*dx 

        self.F += stab
        # self.F_sans_tf += stab

        if self.use_25d_model and self.dom.dim == 2 :
            if self.dom.dim == 3:
                raise ValueError("The 2.5D model requires a 2D simulation.")

            self.fprint("Using 2.5D model")
            dudx = Dx(self.u_next[0], 0)
            dvdy = Dx(self.u_next[1], 1)

            if inflow_angle is None:
                term25 = dvdy*q*dx
            else:
                term25 = (abs(sin(inflow_angle))*dudx*q + abs(cos(inflow_angle))*dvdy*q)*dx

            self.F -= term25

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

    def ComputeFunctional(self,inflow_angle=None):
        self.fprint("Setting Up Taylor-Hood Problem",special="header")

        ### These constants will be moved into the params file ###
        f = Constant((0.0,)*self.dom.dim)
        vonKarman=0.41
        eps=Constant(0.01)
        nu = self.viscosity


        self.fprint("Viscosity:         {:1.2e}".format(float(self.viscosity)))
        self.fprint("Max Mixing Length: {:1.2e}".format(float(self.lmax)))

        ### Create the test/trial/functions ###
        self.up_k = Function(self.fs.W)
        self.u_k,self.p_k = split(self.up_k)
        v,q = TestFunctions(self.fs.W)

        ### Set the initial guess ###
        ### (this will become a separate function.)
        self.up_k.assign(self.bd.u0)

        ### Calculate the stresses and viscosities ###
        S = sqrt(2.*inner(0.5*(grad(self.u_k)+grad(self.u_k).T),0.5*(grad(self.u_k)+grad(self.u_k).T)))

        ### Create l_mix based on distance to the ground ###
        if self.dom.dim == 3:
            ### https://doi.org/10.5194/wes-4-127-2019###
            l_mix = Function(self.fs.Q)
            l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:],(1.+np.divide(vonKarman*self.bd.depth.vector()[:],self.lmax)))
        else:
            l_mix = Constant(vonKarman*self.farm.HH[0]/(1+(vonKarman*self.farm.HH[0]/self.lmax)))

        ### Calculate nu_T
        self.nu_T=l_mix**2.*S

        ### Create the turbine force ###
        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)

        ### Create the functional ###
        self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx + inner(self.tf,v)*dx 

        if self.use_25d_model:
            if self.dom.dim == 3:
                raise ValueError("The 2.5D model requires a 2D simulation.")

            self.fprint("Using 2.5D model")
            dudx = Dx(self.u_k[0], 0)
            dvdy = Dx(self.u_k[1], 1)

            if inflow_angle is None:
                term25 = dvdy*q*dx
            else:
                term25 = (abs(sin(inflow_angle))*dudx*q + abs(cos(inflow_angle))*dvdy*q)*dx

            self.F -= term25

        self.fprint("Taylor-Hood Problem Setup",special="footer")

# ================================================================

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

        ### Create Functional ###
        self.ComputeFunctional()

    def ComputeFunctional(self,inflow_angle=None):
        # ================================================================

        # Define fluid properties
        # FIXME: These should probably be set in params.yaml input filt
        # nu = 1/10000
        rho = 1
        nu_c = Constant(self.viscosity)
        rho_c = Constant(rho)

        # Define time step size (this value is used only for step 1 if adaptive timestepping is used)
        # FIXME: change variable name to avoid confusion within dolfin adjoint
        self.dt = 0.1*self.dom.mesh.hmin()/self.bd.HH_vel
        # self.dt = 0.04
        self.dt_c  = Constant(self.dt)

        self.fprint("Viscosity: {:1.2e}".format(float(self.viscosity)))
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
        self.u_k.assign(self.bd.bc_velocity)
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

        # FIXME: This up_k function is only present to avoid errors  
        # during assignments in GenericSolver.__init__

        # Create the combined function space
        self.up_k = Function(self.fs.W)

        # Create the turbine force
        # FIXME: Should this be set by a numpy array operation or a fenics function?
        # self.tf = self.farm.TurbineForce(self.fs, self.dom.mesh, self.u_k2)
        # self.tf = Function(self.fs.V)

        self.tf = self.ComputeTurbineForce(self.u_k,inflow_angle)
        self.u_k.assign(self.tf)
        self.u_k.assign(self.bd.bc_velocity)


        # self.u_k2.vector()[:] = 0.0
        # self.u_k1.vector()[:] = 0.0

        # ================================================================

        # Define variational problem for step 1: tentative velocity
        # F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
        #    + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
        #    + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
        #    + dot(nabla_grad(self.p_k1), v)*dx \
        #    - dot(-self.tf, v)*dx

        F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
           + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
           + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
           + dot(nabla_grad(self.p_k1), v)*dx \
           - dot(self.tf, v)*dx

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

# # ================================================================

#     def UpdateActuatorLineControls(self, c_lift = None, c_drag = None):

#         if c_lift is not None:
#             cl = np.array(c_lift, dtype = float)
#         if c_drag is not None:
#             cd = np.array(c_drag, dtype = float)

#         for k in range(self.num_blade_segments):
#             self.mcl[k] = Constant(cl[k])
#             self.mcd[k] = Constant(cd[k])

# # ================================================================

