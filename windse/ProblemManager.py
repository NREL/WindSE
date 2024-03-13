"""
The ProblemManager contains all of the 
different classes of problems that windse can solve
"""

import __main__
import os
from pyadjoint.tape import no_annotations

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    from dolfin import *
    import numpy as np
    import time
    import scipy.interpolate as interp
    import glob

    ### Import the cumulative parameters ###
    from windse import windse_parameters
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
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode

        ### Update attributes based on params file ###
        for key, value in self.params["problem"].items():
            setattr(self,key,value)

        self.record_time = self.params["optimization"].get("record_time",0.0)
        if isinstance(self.record_time,str):
            self.record_time = 0.0

        self.extra_kwarg = {}            
        if self.params.dolfin_adjoint:
            self.extra_kwarg["annotate"] = False

        # Setup body force
        self.mbody_force = Constant(self.body_force)


    @no_annotations
    def DebugOutput(self):
        if self.debug_mode:
            # integral of nu_t
            int_nut = assemble(self.nu_T*dx)/self.dom.volume
            self.tag_output("int_nu_T", int_nut)

            if self.dom.dim == 3:
                e1 = Constant((1,0,0)); e2 = Constant((0,1,0)); e3 = Constant((0,0,1));
            else:
                e1 = Constant((1,0)); e2 = Constant((0,1));

            if self.farm.use_local_tf_dx:
                int_tf_x = 0
                int_tf_y = 0
                int_tf_z = 0
                for i in range(len(self.farm.tf_list)):
                    tf = self.farm.tf_list[i]
                    local_dx = self.farm.local_dx(i+1)

                    int_tf_x += assemble(inner(tf,e1)*local_dx)/self.dom.volume
                    int_tf_y += assemble(inner(tf,e2)*local_dx)/self.dom.volume
                    if self.dom.dim == 3:
                        int_tf_z += assemble(inner(tf,e3)*local_dx)/self.dom.volume
            else:
                tf = 0
                for i in range(len(self.farm.tf_list)):
                    tf += self.farm.tf_list[i]

                int_tf_x = assemble(inner(tf,e1)*dx)/self.dom.volume
                int_tf_y = assemble(inner(tf,e2)*dx)/self.dom.volume
                if self.dom.dim == 3:
                    int_tf_z = assemble(inner(tf,e3)*dx)/self.dom.volume

            self.tag_output("int_tf_x", int_tf_x)
            self.tag_output("int_tf_y", int_tf_y)

            if self.dom.dim == 3:
                self.tag_output("int_tf_z", int_tf_z)

                if self.farm.turbine_type == 'line':
                    chord = []
                    cl = []
                    cd = []
                    num_actuator_nodes = []
                    for i in range(self.farm.numturbs):
                        chord.append(self.farm.turbines[i].chord)
                        cl.append(self.farm.turbines[i].cl)
                        cd.append(self.farm.turbines[i].cd)
                        num_actuator_nodes.append(self.farm.turbines[i].num_actuator_nodes)

                    self.tag_output("min_chord", np.min(chord))
                    self.tag_output("max_chord", np.max(chord))
                    self.tag_output("avg_chord", np.mean(chord))
                    self.tag_output("min_cl", np.min(cl))
                    self.tag_output("max_cl", np.max(cl))
                    self.tag_output("avg_cl", np.mean(cl))
                    self.tag_output("min_cd", np.min(cd))
                    self.tag_output("max_cd", np.max(cd))
                    self.tag_output("avg_cd", np.mean(cd))
                    self.tag_output("num_actuator_nodes", np.mean(num_actuator_nodes))
                

    def ComputeTurbineForceTerm(self,u,v,inflow_angle,**kwargs):
        tf_start = time.time()
        self.fprint("Calculating Turbine Force",special="header")

        # Compute the relative yaw angle
        if inflow_angle is None:
            inflow_angle = self.dom.inflow_angle

        self.fprint('Computing turbine forces using %s' % (self.farm.turbine_type.upper()))

        # Create a list of turbine force function and domain of integration for each turbine
        if self.farm.turbine_type == "disabled" or self.farm.numturbs == 0:

            # if there are no turbine return an zero force term
            tf_term = Function(self.fs.V)*dx
        else:

            # compute tf and dx for each turbine
            tf_term = self.farm.compute_turbine_force(u,v,inflow_angle,self.fs,**kwargs)

        tf_stop = time.time()
        self.fprint("Turbine Force Calculated: {:1.2f} s".format(tf_stop-tf_start),special="footer")
        return tf_term

    def ComputeTurbulenceModel(self, u):
        self.fprint(f"Using Turbulence Model: {self.turbulence_model}")
        if self.turbulence_model is not None:
            # Calculate eddy viscosity, if using
            # self.turbulence_model = 'smagorinsky'
            # self.turbulence_model = 'mixing_length'
            # self.turbulence_model = None

            # Strain rate tensor, 0.5*(du_i/dx_j + du_j/dx_i)
            Sij = sym(grad(u))

            # sqrt(Sij*Sij)
            strainMag = (2.0*inner(Sij, Sij))**0.5

            if self.turbulence_model == 'smagorinsky':
                # Smagorinsky constant, typically around 0.17
                Cs = 0.17

                # Define filter scale (all models use this)
                filter_scale = CellVolume(self.dom.mesh)**(1.0/self.dom.dim)

                # Eddy viscosity
                nu_T = Cs**2 * filter_scale**2 * strainMag

            elif self.turbulence_model == 'mixing_length':
                # von Karman constant
                vonKarman = 0.41

                if self.dom.dim == 3:
                    depth = self.bd.depth.vector()[:]
                else:
                    if self.farm.numturbs == 0:
                        depth = self.params["boundary_conditions"]["vel_height"]
                    else:
                        depth = self.farm.turbines[0].HH

                numer = vonKarman*depth
                denom = 1.0 + vonKarman*depth/self.lmax

                l_mix = Function(self.fs.Q)
                l_mix.vector()[:] = numer/denom
                l_mix.rename("l_mix","l_mix")

                # Eddy viscosity
                nu_T = l_mix**2 * strainMag

        else:
            nu_T = Constant(0)

        return nu_T


    ############################################
    ############################################
    ### this is a hack to accommodate the fact that
    ### alm values are not stored in the wind_farm
    ### object. eventually we might want to move
    ### it there.
    ############################################
    ############################################
    def CopyALMtoWindFarm(self):
        self.farm.mcl = self.mcl
        self.farm.mcd = self.mcd
        self.farm.mchord = self.mchord
        self.farm.cl = self.cl
        self.farm.cd = self.cd
        self.farm.chord = self.chord
        self.farm.num_actuator_nodes = self.num_actuator_nodes

    def SimpleControlUpdate(self):
        self.u_k, self.p_k = split(self.up_k)
        self.farm.SimpleControlUpdate()


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

    def UpdateActuatorLineControls(self, c_lift = None, c_drag = None, chord = None, yaw = None, turb_i = 0):

        if c_lift is not None:
            cl = np.array(c_lift, dtype = float)
            self.cl[turb_i] = cl
            for k in range(self.num_actuator_nodes):
                self.mcl[turb_i][k] = Constant(cl[k])
        if c_drag is not None:
            cd = np.array(c_drag, dtype = float)
            self.cd[turb_i] = cd
            for k in range(self.num_actuator_nodes):
                self.mcd[turb_i][k] = Constant(cd[k])
        if chord is not None:
            chord = np.array(chord, dtype = float)
            self.chord[turb_i] = chord
            for k in range(self.num_actuator_nodes):
                self.mchord[turb_i][k] = Constant(chord[k])
        if yaw is not None:
            yaw = float(yaw)
            self.farm.yaw[turb_i] = yaw
            self.farm.myaw[turb_i] = Constant(yaw)
        

        self.CopyALMtoWindFarm()


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
        self.ComputeFunctional(self.dom.inflow_angle)
        self.DebugOutput()


    def ComputeFunctional(self,inflow_angle):
        self.fprint("Setting Up Stabilized Problem",special="header")

        ### Create the test/trial/functions ###
        self.up_k = Function(self.fs.W)
        self.u_k,self.p_k = split(self.up_k)
        v,q = TestFunctions(self.fs.W)

        ### Set the x scaling ###
        Sx = self.dom.xscale

        ### Set the initial guess ###
        u0 = self.bd.u0
        self.up_k.assign(u0)

        # mem0=memory_usage()[0]
        # mem_out, self.tf = memory_usage((self.ComputeTurbineForce,(self.u_k,inflow_angle),{}),max_usage=True,retval=True,max_iterations=1)
        # self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))
        tf_term = self.ComputeTurbineForceTerm(self.u_k,v,inflow_angle)

        # set_log_level(LogLevel.DEBUG)
        # tic = time.time()
        # assemble(turbine_force_term)
        # toc = time.time()
        # print(f"assemble time: {toc-tic} s")
        # exit()


        ### These constants will be moved into the params file ###
        f = Constant((0.0,)*self.dom.dim)
        f.rename("f","f")
        
        nu = self.viscosity
        vonKarman=0.41
        eps=Constant(self.stability_eps)
        eps.rename("eps","eps")

        self.fprint("Viscosity:                 {:1.2e}".format(float(self.viscosity)))
        self.fprint("Max Mixing Length:         {:1.2e}".format(float(self.lmax)))
        self.fprint("Stabilization Coefficient: {:1.2e}".format(float(eps)))

        ### Calculate nu_T
        self.nu_T=self.ComputeTurbulenceModel(self.u_k)
        # self.nu_T=Constant(0.0)
        self.ReyStress=self.nu_T*grad(self.u_k)
        self.vertKE= self.ReyStress[0,2]*self.u_k[0]

        ### Create the functional ###
        self.F = inner(grad(self.u_k)*self.u_k, v)*dx
        self.F +=   Sx*Sx*(nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx
        self.F += - inner(div(v),self.p_k)*dx
        self.F += - inner(div(self.u_k),q)*dx
        self.F += - inner(f,v)*dx
        self.F += - tf_term 


        # Add body force to functional
        if abs(float(self.mbody_force)) >= 1e-14:
            self.fprint("Using Body Force")
            self.F += inner(-self.mbody_force*self.bd.inflow_unit_vector,v)*dx

        ################ THIS IS A CHEAT ####################\
        if self.use_corrective_force:
            self.fprint("Using Corrective Force")
            extra_S = sqrt(2*inner(0.5*(grad(self.bd.bc_velocity)+grad(self.bd.bc_velocity).T),0.5*(grad(self.bd.bc_velocity)+grad(self.bd.bc_velocity).T)))
            extra_l_mix = Function(self.fs.Q)
            extra_l_mix.vector()[:] = np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,(1.+np.divide(vonKarman*self.bd.depth.vector()[:]/Sx,self.lmax)))
            extra_nu_T = extra_l_mix**2.*extra_S
            extra_DP =dot(self.bd.u0,grad(self.bd.u0)) - div((nu+extra_nu_T)*grad(self.bd.bc_velocity))

            self.F += inner(extra_DP,v)*dx
        ########################################################

        # self.F_sans_tf =  (1.0)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx
        # self.F = inner(grad(self.u_k)*self.u_k, v)*dx + (nu+self.nu_T)*inner(grad(self.u_k), grad(v))*dx - inner(div(v),self.p_k)*dx - inner(div(self.u_k),q)*dx - inner(f,v)*dx - inner(self.tf*(self.u_k[0]**2+self.u_k[1]**2),v)*dx 
        # stab_sans_tf = - eps*inner(grad(q), grad(self.p_k))*dx 
        # self.F_sans_tf += stab

        ### Add in the Stabilizing term ###
        if abs(float(eps)) >= 1e-14:
            self.fprint("Using Stabilization Term")
            stab = - eps*inner(grad(q), grad(self.p_k))*dx - eps*inner(grad(q), dot(grad(self.u_k), self.u_k))*dx 
            self.F += stab


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


        self.first_loop = True

        ### Create Functional ###
        self.ComputeFunctional(self.dom.inflow_angle)
        self.DebugOutput()

    def ComputeFunctional(self,inflow_angle):
        self.fprint("Setting Up Taylor-Hood Problem",special="header")

        ### These constants will be moved into the params file ###
        f = Constant((0.0,)*self.dom.dim)
        vonKarman=0.41
        eps=Constant(0.000001)
        nu = self.viscosity


        self.fprint("Viscosity:         {:1.2e}".format(float(self.viscosity)))
        self.fprint("Max Mixing Length: {:1.2e}".format(float(self.lmax)))


        ### Create the test/trial/functions ###
        self.v,self.q = TestFunctions(self.fs.W)

        ### Set the initial guess ###
        ### (this will become a separate function.)
        # if self.first_loop:
        #     self.up_k = Function(self.fs.W)
        #     self.up_k.assign(self.bd.u0)
        #     self.first_loop = False
        # else:
        #     temp_up = self.up_k.copy()
        self.up_k = Function(self.fs.W)
        self.up_k.assign(self.bd.u0)
        self.u_k,self.p_k = split(self.up_k)

        ### Calculate nu_T
        self.nu_T=self.ComputeTurbulenceModel(self.u_k)

        ### Create the turbine force ###
        tf_term = self.ComputeTurbineForceTerm(self.u_k,self.v,inflow_angle)

        ### Create the functional ###
        self.F = inner(grad(self.u_k)*self.u_k, self.v)*dx
        self.F +=   (nu+self.nu_T)*inner(grad(self.u_k), grad(self.v))*dx
        self.F += - inner(div(self.v),self.p_k)*dx
        self.F += - inner(div(self.u_k),self.q)*dx
        # self.F += - inner(f,v)*dx
        self.F += - tf_term 

        # ### Add in the Stabilizing term ###
        # if abs(float(eps)) >= 1e-14:
        #     self.fprint("Using Stabilization Term")
        #     stab = - eps*inner(grad(q), grad(self.p_k))*dx - eps*inner(grad(q), dot(grad(self.u_k), self.u_k))*dx 
        #     self.F += stab

        # Add body force to functional
        if abs(float(self.mbody_force)) >= 1e-14:
            self.fprint("Using Body Force")
            self.F += inner(-self.mbody_force*self.bd.inflow_unit_vector,self.v)*dx

        if self.use_25d_model:
            if self.dom.dim == 3:
                raise ValueError("The 2.5D model requires a 2D simulation.")

            self.fprint("Using 2.5D model")
            dudx = Dx(self.u_k[0], 0)
            dvdy = Dx(self.u_k[1], 1)

            if inflow_angle is None:
                term25 = dvdy*self.q*dx
            else:
                term25 = (abs(sin(inflow_angle))*dudx*self.q + abs(cos(inflow_angle))*dvdy*self.q)*dx

            self.F -= term25

        self.fprint("Taylor-Hood Problem Setup",special="footer")

# ================================================================

class IterativeSteady(GenericProblem):
    """
    The IterativeSteady sets up everything required for solving Navier-Stokes using
    the SIMPLE algorithm

    Args: 
        domain (:meth:`windse.DomainManager.GenericDomain`): a windse domain object.
        windfarm (:meth:`windse.WindFarmManager.GenericWindFarmm`): a windse windfarm object.
        function_space (:meth:`windse.FunctionSpaceManager.GenericFunctionSpace`): a windse function space object.
        boundary_conditions (:meth:`windse.BoundaryManager.GenericBoundary`): a windse boundary object.
    """
    def __init__(self, domain, windfarm, function_space, boundary_conditions):
        super(IterativeSteady, self).__init__(domain, windfarm, function_space, boundary_conditions)
        self.fprint("Setting Up *Iterative* Steady Problem", special="header")

        ### Create Functional ###
        self.ComputeFunctional(self.dom.inflow_angle)

    def ComputeFunctional(self,inflow_angle):
        # ================================================================

        # Define fluid properties
        # FIXME: These should probably be set in params.yaml input filt
        # nu = 1/10000
        nu = Constant(self.viscosity)

        # Trial functions for velocity and pressure
        u = TrialFunction(self.fs.V)
        p = TrialFunction(self.fs.Q)

        # Test functions for velocity and pressure
        v = TestFunction(self.fs.V)
        q = TestFunction(self.fs.Q)

        ### Define Velocity Functions ###
        self.u_k = Function(self.fs.V, name="u_k")
        self.u_hat = Function(self.fs.V)
        self.u_s = Function(self.fs.V)
        self.du = Function(self.fs.V)
        self.u_k_old = Function(self.fs.V)

        ### Define Pressure Functions ###
        self.p_k = Function(self.fs.Q, name="p_k")
        self.p_s = Function(self.fs.Q)
        self.dp = Function(self.fs.Q)
        self.p_k_old = Function(self.fs.Q)

        U_CN = 0.5*(u + self.u_k) 

        # Adams-Bashforth velocity
        # U_AB = 1.5*u_k - 0.5*u_k_old # Time level k+1/2
        U_AB = 2.0*self.u_k - 1.0*self.u_k_old # Time level k+1

        # Compute eddy viscosity
        self.nu_T=self.ComputeTurbulenceModel(U_AB)

        # ================================================================

        # FIXME: This up_k function is only present to avoid errors  
        # during assignments in GenericSolver.__init__

        # Create the combined function space
        self.up_k = Function(self.fs.W)

        # Create the turbine force
        # FIXME: Should this be set by a numpy array operation or a fenics function?
        # self.tf = self.farm.TurbineForce(self.fs, self.dom.mesh, self.u_k2)
        # self.tf = Function(self.fs.V)

        tf_term = self.ComputeTurbineForceTerm(self.u_k,v,inflow_angle)
        # self.u_k.assign(self.bd.bc_velocity)

        # self.u_k2.vector()[:] = 0.0
        # self.u_k1.vector()[:] = 0.0

        # ================================================================

        # Solve for u_hat, a velocity estimate which doesn't include pressure gradient effects
        F1 = inner(dot(self.u_k, nabla_grad(u)), v)*dx \
           + (nu+self.nu_T)*inner(grad(u), grad(v))*dx \
           - tf_term 

        # Add body force to functional
        F1 += inner(-self.mbody_force*self.bd.inflow_unit_vector,v)*dx

        self.F1_lhs = lhs(F1)
        self.F1_rhs = rhs(F1)


        # Use u_hat to solve for the pressure field
        self.F2_lhs = inner(grad(p), grad(q))*dx
        self.F2_rhs = - div(self.u_hat)*q*dx

        self.dt_1 = Constant(1) # Becomes (1/dt, 0) for (unsteady, steady) state
        self.dt_2 = Constant(1) # Becomes (1/dt, 1) for (unsteady, steady) state
        self.dt_3 = Constant(1) # Becomes (  dt, 1) for (unsteady, steady) state

        # Solve for u_star, a predicted velocity which includes the pressure gradient
        F3 = inner(dot(self.u_k, nabla_grad(u)), v)*dx \
           + (nu+self.nu_T)*inner(grad(u), grad(v))*dx \
           - tf_term \
           + inner(grad(self.p_k), v)*dx \
           + self.dt_1*inner(u - self.u_k, v)*dx

        # Add body force to functional
        F3 += inner(-self.mbody_force*self.bd.inflow_unit_vector,v)*dx

        self.F3_lhs = lhs(F3)
        self.F3_rhs = rhs(F3)


        # Define variational problem for step 2: pressure correction
        self.F4_lhs = inner(grad(p), grad(q))*dx
        self.F4_rhs = - self.dt_2*div(self.u_s)*q*dx


        # Define variational problem for step 3: velocity update
        self.F5_lhs = inner(u, v)*dx
        self.F5_rhs = - self.dt_3*inner(grad(self.dp), v)*dx

        # ================================================================

        self.fprint("*Iterative* Steady Problem Setup",special="footer")

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
        self.ComputeFunctional(self.dom.inflow_angle)
        self.DebugOutput()

    def ComputeFunctional(self,inflow_angle):
        # ================================================================

        # Define fluid properties
        # FIXME: These should probably be set in params.yaml input filt
        # nu = 1/10000
        rho = 1
        nu_c = Constant(self.viscosity, name="viscosity")
        rho_c = Constant(rho, name="rho")

        # Define time step size (this value is used only for step 1 if adaptive timestepping is used)
        # FIXME: change variable name to avoid confusion within dolfin adjoint
        self.dt = 0.1*self.dom.global_hmin/self.bd.HH_vel
        # self.dt = 0.05
        self.dt_c  = Constant(self.dt, name="dt_c")

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
        self.u_k = Function(self.fs.V, name="u_k")
        self.u_k1 = Function(self.fs.V, name="u_k1")
        self.u_k2 = Function(self.fs.V, name="u_k2")

        # Seed previous velocity fields with the chosen initial condition
        self.u_k.assign(self.bd.bc_velocity)
        self.u_k1.assign(self.bd.bc_velocity)
        self.u_k2.assign(self.bd.bc_velocity)

        # Calculate Reynolds stress 
        self.uk_sum = Function(self.fs.V, name="uk_sum")
        self.uk_sum.assign(self.dt_c*self.u_k)
        self.vertKE = Function(self.fs.Q, name="vertKE")

        # Define functions for pressure solutions
        # >> _k = current (step k)
        # >> _k1 = previous (step k-1)
        self.p_k  = Function(self.fs.Q, name="p_k")
        self.p_k1 = Function(self.fs.Q, name="p_k1")

        # Seed previous pressure fields with the chosen initial condition
        self.p_k1.assign(self.bd.bc_pressure)

        # ================================================================

        # Crank-Nicolson velocity
        U_CN  = 0.5*(u + self.u_k1)

        # Adams-Bashforth projected velocity
        U_AB = 1.5*self.u_k1 - 0.5*self.u_k2

        # ================================================================
        self.nu_T = self.ComputeTurbulenceModel(U_AB)

        # ================================================================

        # FIXME: This up_k function is only present to avoid errors  
        # during assignments in GenericSolver.__init__

        # Create the combined function space
        self.up_k = Function(self.fs.W, name="up_k")

        # Create the turbine force
        # FIXME: Should this be set by a numpy array operation or a fenics function?
        # self.tf = self.farm.TurbineForce(self.fs, self.dom.mesh, self.u_k2)
        # self.tf = Function(self.fs.V)

        tf_term = self.ComputeTurbineForceTerm(self.u_k,v,inflow_angle,simTime=0.0,simTime_prev=None, dt=self.dt)

        self.u_k.assign(self.bd.bc_velocity)

        # Only the actuator lines point "upstream" against the flow
        # all other actuator forces need to be reversed to follow
        # the convention used in the unsteady problem
        # FIXME: We should be consistent about which way the turbine
        # forces are oriented.
        # if self.farm.turbine_method != "alm":
        #     self.tf *= -1.0

        # self.u_k2.vector()[:] = 0.0
        # self.u_k1.vector()[:] = 0.0

        # ================================================================

        # Define variational problem for step 1: tentative velocity
        # F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
        #    + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
        #    + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
        #    + dot(nabla_grad(self.p_k1), v)*dx \
        #    - dot(-self.tf, v)*dx

        # F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
        #    + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
        #    + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
        #    + dot(nabla_grad(self.p_k1), v)*dx \
        #    - dot(self.tf, v)*dx

        F1 = (1.0/self.dt_c)*inner(u - self.u_k1, v)*dx \
           + inner(dot(U_AB, nabla_grad(U_CN)), v)*dx \
           + (nu_c+self.nu_T)*inner(grad(U_CN), grad(v))*dx \
           + inner(grad(self.p_k1), v)*dx \
           - tf_term

        self.a1 = lhs(F1)
        self.L1 = rhs(F1)

        # Define variational problem for step 2: pressure correction
        # self.a2 = dot(nabla_grad(p), nabla_grad(q))*dx
        # self.L2 = dot(nabla_grad(self.p_k1), nabla_grad(q))*dx - (1.0/self.dt_c)*div(self.u_k)*q*dx
        self.a2 = inner(grad(p), grad(q))*dx
        self.L2 = inner(grad(self.p_k1), grad(q))*dx - (1.0/self.dt_c)*div(self.u_k)*q*dx

        # phi = p - self.p_k
        # F2 = inner(grad(q), grad(phi))*dx - (1.0/self.dt_c)*div(u_k)*q*dx
        # self.a2 = lhs(F2)
        # self.L2 = rhs(F2)

        # Define variational problem for step 3: velocity update
        # self.a3 = dot(u, v)*dx
        # self.L3 = dot(self.u_k, v)*dx - self.dt_c*dot(nabla_grad(self.p_k - self.p_k1), v)*dx
        self.a3 = inner(u, v)*dx
        self.L3 = inner(self.u_k, v)*dx - self.dt_c*inner(grad(self.p_k - self.p_k1), v)*dx

        # F3 = inner(u, v)*dx - inner(self.u_k, v)*dx + self.dt_c*inner(phi, v)*dx
        # self.a3 = lhs(F3)
        # self.L3 = rhs(F3)


        # ================================================================

        self.fprint("Unsteady Problem Setup",special="footer")

# # ================================================================

#     def UpdateActuatorLineControls(self, c_lift = None, c_drag = None):

#         if c_lift is not None:
#             cl = np.array(c_lift, dtype = float)
#         if c_drag is not None:
#             cd = np.array(c_drag, dtype = float)

#         for k in range(self.num_actuator_nodes):
#             self.mcl[k] = Constant(cl[k])
#             self.mcd[k] = Constant(cd[k])

# # ================================================================

