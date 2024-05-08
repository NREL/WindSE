"""
The SolverManager contains all the different ways to solve problems generated
in windse
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
    from sys import platform
    import time
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
    # from memory_profiler import memory_usage
    from scipy.special import gamma
    from mpmath import hyper
    import cProfile
    import pstats

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Import some helper functions
    from windse.helper_functions import meteor_to_math, math_to_meteor

    ### Check if we need dolfin_adjoint ###
    if windse_parameters.dolfin_adjoint:
        from dolfin_adjoint import *
        from windse.blocks import blockify, MarkerBlock, ControlUpdaterBlock
        from pyadjoint import AdjFloat

    ### Import objective functions ###
    import windse.objective_functions as obj_funcs

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt

    ### Check to see if the latest compiler representation is installed ###
    try:
        import tsfc
        default_representation = 'tsfc'
    except ImportError as e:
        default_representation = 'uflacs'

    ### Improved dolfin parameters ###
    parameters["std_out_all_processes"] = False;
    parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -fno-math-errno -march=native'        
    parameters["form_compiler"]["optimize"]     = True
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters['form_compiler']['representation'] = default_representation
    parameters['form_compiler']['quadrature_degree'] = windse_parameters["function_space"]["quadrature_degree"]
    
class GenericSolver(object):
    """
    A GenericSolver contains on the basic functions required by all solver objects.
    """
    def __init__(self,problem):
        self.params = windse_parameters
        self.problem  = problem
        self.nu_T = self.problem.nu_T
        self.first_save = True
        self.fprint = self.params.fprint
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode
        self.simTime = 0.0
        self.simTime_prev = None
        self.iter_val = 0.0
        self.pow_saved = False

        ### Update attributes based on params file ###
        for key, value in self.params["solver"].items():
            setattr(self,key,value)

        ### Update attributes based on params file ###
        for key, value in self.params["optimization"].items():
            setattr(self,key,value)

        self.extra_kwarg = {}            
        if self.params.dolfin_adjoint:
            self.extra_kwarg["annotate"] = False

        J_temp = 0.0
        if self.params["general"]["dolfin_adjoint"]:
            J_temp = AdjFloat(0.0)

        self.optimizing = False
        if self.params.performing_opt_calc or self.save_objective:
            self.optimizing = True
            self.J = J_temp
            self.adj_file = XDMFFile(self.params.folder+"timeSeries/local_adjoint.xdmf")
            self.adj_time_iter = 1
            self.adj_time_list = [0.0]
            self.J_saved = False

            if isinstance(self.opt_turb_id,int):
                self.opt_turb_id = [self.opt_turb_id]
            elif self.opt_turb_id == "all":
                self.opt_turb_id = range(self.problem.farm.numturbs) 
            elif isinstance(self.opt_turb_id, str):
                self.opt_turb_id = [int(self.opt_turb_id)]

        # blockify custom functions so dolfin adjoint can track them
        if self.params.performing_opt_calc:
        #     self.marker = blockify(self.marker,MarkerBlock)
            self.control_updater = blockify(self.control_updater,ControlUpdaterBlock)

        #Check if we need to save the power output
        if self.save_power:
            self.J = J_temp
            self.J_saved = False

    @no_annotations
    def DebugOutput(self,t=None,i=None):
        if self.debug_mode:
            if self.problem.dom.dim == 3:
                ux, uy, uz = self.problem.u_k.split(True)
            else:
                ux, uy = self.problem.u_k.split(True)

            if t is None:
                suffix = ""
            else:
                suffix = "_"+repr(i)
                self.tag_output("time"+suffix,t)

            self.tag_output("min_x_vel"+suffix,ux.vector().min()) # probably not the fastest way to get the average velocity
            self.tag_output("max_x_vel"+suffix,ux.vector().max()) # probably not the fastest way to get the average velocity
            self.tag_output("avg_x_vel"+suffix,assemble(ux*dx)/self.problem.dom.volume) # probably not the fastest way to get the average velocity
            self.tag_output("min_y_vel"+suffix,uy.vector().min()) # probably not the fastest way to get the average velocity
            self.tag_output("max_y_vel"+suffix,uy.vector().max()) # probably not the fastest way to get the average velocity
            self.tag_output("avg_y_vel"+suffix,assemble(uy*dx)/self.problem.dom.volume) # probably not the fastest way to get the average velocity
            if self.problem.dom.dim == 3:
                self.tag_output("min_z_vel"+suffix,uz.vector().min()) # probably not the fastest way to get the average velocity
                self.tag_output("max_z_vel"+suffix,uz.vector().max()) # probably not the fastest way to get the average velocity
                self.tag_output("avg_z_vel"+suffix,assemble(uy*dx)/self.problem.dom.volume) # probably not the fastest way to get the average velocity



    @no_annotations
    def Save(self,val=0):
        """
        This function saves the mesh and boundary markers to output/.../solutions/
        """
        u = self.problem.u_k
        p = self.problem.p_k

        u.vector()[:]=u.vector()[:]/self.problem.dom.xscale
        self.problem.dom.mesh.coordinates()[:]=self.problem.dom.mesh.coordinates()[:]/self.problem.dom.xscale

        if self.first_save:
            self.u_file = self.params.Save(u,"velocity",subfolder="solutions/",val=val)
            self.p_file = self.params.Save(p,"pressure",subfolder="solutions/",val=val)
            self.nuT_file = self.params.Save(self.nu_T,"eddy_viscosity",subfolder="solutions/",val=val)
            if self.problem.dom.dim == 3:
                self.ReyStress_file = self.params.Save(self.ReyStress,"Reynolds_stresses",subfolder="solutions/",val=val)
                self.vertKE_file = self.params.Save(self.vertKE,"Vertical KE",subfolder="solutions/",val=val)
            self.first_save = False
        else:
            self.params.Save(u,"velocity",subfolder="solutions/",val=val,file=self.u_file)
            self.params.Save(p,"pressure",subfolder="solutions/",val=val,file=self.p_file)
            self.params.Save(self.nu_T,"eddy_viscosity",subfolder="solutions/",val=val,file=self.nuT_file)
            if self.problem.dom.dim == 3:
                self.params.Save(self.ReyStress,"Reynolds_stresses",subfolder="solutions/",val=val,file=self.ReyStress_file)
                self.params.Save(self.vertKE,"Vertical KE",subfolder="solutions/",val=val,file=self.vertKE_file)
        u.vector()[:]=u.vector()[:]*self.problem.dom.xscale
        self.problem.dom.mesh.coordinates()[:]=self.problem.dom.mesh.coordinates()[:]*self.problem.dom.xscale

    def ChangeWindSpeed(self,inflow_speed):
        """
        This function recomputes all necessary components for a new wind direction

        Args: 
            theta (float): The new wind angle in radians
        """
        self.problem.ChangeWindSpeed(inflow_speed)

    def ChangeWindAngle(self,inflow_angle):
        """
        This function recomputes all necessary components for a new wind direction

        Args: 
            theta (float): The new wind angle in radians
        """
        self.problem.ChangeWindAngle(inflow_angle)

    # def SavePower(self,inflow_angle=0.0):

    #     J_list=np.zeros(self.problem.farm.numturbs+1)

    #     if self.problem.farm.actuator_disks_list is not None:
    #         for i in range(self.problem.farm.numturbs):
    #             yaw = self.problem.farm.myaw[i]+inflow_angle
    #             tf1 = self.problem.farm.actuator_disks_list[i] * cos(yaw)**2
    #             tf2 = self.problem.farm.actuator_disks_list[i] * sin(yaw)**2
    #             tf3 = self.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
    #             tf = tf1*self.u_k[0]**2+tf2*self.u_k[1]**2+tf3*self.u_k[0]*self.u_k[1]
    #             J_list[i] = assemble(dot(tf,self.u_k)*dx,**self.extra_kwarg)
    #     J_list[-1]=sum(J_list)

    #     folder_string = self.params.folder+"/data/"
    #     if not os.path.exists(folder_string): os.makedirs(folder_string)

    #     if self.J_saved:
    #         f = open(folder_string+"power_data.txt",'ab')
    #     else:
    #         f = open(folder_string+"power_data.txt",'wb')
    #         self.J_saved = True

    #     np.savetxt(f,[J_list])
    #     f.close()

    @no_annotations
    def EvaulatePowerFunctional(self):
        kwargs = {
            "iter_val": self.iter_val, 
            "simTime": self.simTime
        }
        out = self.problem.farm.save_power(self.problem.u_k,self.problem.dom.inflow_angle, **kwargs)
        return out


    def EvaluateObjective(self,output_name="objective_data",opt_iter=-1):
        self.fprint("Evaluating Objective Data",special="header")
        start = time.time()

        first_call = True
        if self.J_saved:
            first_call = False

        annotate = self.params.dolfin_adjoint 

        ### Iterate over objectives ###
        obj_list = [opt_iter, self.iter_val, self.simTime]
        for key, obj_kwargs in self.objective_type.items():
            if isinstance(key, int):
                objective_name = obj_kwargs["type"]
            else:
                objective_name = key
            objective_func = obj_funcs.objective_functions[objective_name]
            args = (self, (self.problem.dom.inflow_angle))
            kwargs = {"first_call": first_call, "annotate": annotate}
            kwargs.update(obj_kwargs)
            out = obj_funcs._annotated_objective(objective_func, *args, **kwargs)
            obj_list.append(out)
        J = obj_list[3] #grab first objective 

        # ### Flip the sign because the objective is minimized but these values are maximized
        # for i in range(1,len(obj_list)):
        #     obj_list[i] = -obj_list[i]

        ### Save to csv ###
        if self.J_saved:
            self.params.save_csv(output_name,data=[obj_list],subfolder=self.params.folder+"data/",mode='a')
        else:
            ### Generate the header ###
            header = "Opt_iter, Iter_Val, Time, "
            for key, val in self.objective_type.items():
                if isinstance(key, int):
                    name = val["type"]
                else:
                    name = key
                header += name + ", "
            header = header[:-2]

            self.params.save_csv(output_name,header=header,data=[obj_list],subfolder=self.params.folder+"data/",mode='w')
            self.J_saved = True

        stop = time.time()
        self.fprint("Complete: {:1.2f} s".format(stop-start),special="footer")
        return J

    def marker(self, u, simTime, adj_tape_file):
        return u

    def control_updater(self, J, problem, time=None):
        return J





class SteadySolver(GenericSolver):
    """
    This solver is for solving the steady state problem

    Args: 
        problem (:meth:`windse.ProblemManager.GenericProblem`): a windse problem object.
    """
    def __init__(self,problem):
        super(SteadySolver, self).__init__(problem)
        self.u_k,self.p_k = split(self.problem.up_k)

    def Solve(self):
        """
        This solves the problem setup by the problem object.
        """

        ### Save Files before solve ###
        self.fprint("Saving Input Data",special="header")
        if "mesh" in self.params.output:
            self.problem.dom.Save(val=self.iter_val)
        if "initial_guess" in self.params.output:
            self.problem.bd.SaveInitialGuess(val=self.iter_val)
        if "height" in self.params.output and self.problem.dom.dim == 3:
            self.problem.bd.SaveHeight(val=self.iter_val)
        if "turbine_force" in self.params.output:
            self.problem.farm.save_functions(val=self.iter_val)
        self.fprint("Finished",special="footer")

        # exit()

        ####################################################################
        ### This is the better way to define a nonlinear problem but it 
        ### doesn't play nice with dolfin_adjoint
        ### Define Jacobian ###
        # dU = TrialFunction(self.problem.fs.W)
        # J  = derivative(self.problem.F,  self.problem.up_k, dU)

        # ### Setup nonlinear solver ###
        # nonlinear_problem = NonlinearVariationalProblem(self.problem.F, self.problem.up_k, self.problem.bd.bcs, J)
        # nonlinear_solver  = NonlinearVariationalSolver(nonlinear_problem)

        # ### Set some parameters ###
        # solver_parameters = nonlinear_solver.parameters['newton_solver']['linear_solver'] = 'gmres'

        # # nonlinear_solver.ksp().setGMRESRestart(40)

        # def print_nested_dict(d, indent):
        #     for key, value in d.items():
        #         if hasattr(value, 'items'):
        #             for i in range(indent):
        #                 print('\t', end = '')
        #             print(key, ":")
        #             indent += 1
        #             print_nested_dict(value, indent)
        #             indent -= 1
        #         else:
        #             for i in range(indent):
        #                 print('\t', end = '')
        #             print(key, ":", value)

        # nonlinear_solver.parameters
        # print_nested_dict(nonlinear_solver.parameters, 0)
        # exit()



        # print(type(solver_parameters))

        # info(solver_parameters)
        # solver_parameters["nonlinear_solver"] = "snes"
        # solver_parameters["snes_solver"]["linear_solver"] = "mumps"
        # solver_parameters["snes_solver"]["maximum_iterations"] = 50
        # solver_parameters["snes_solver"]["error_on_nonconvergence"] = False
        # solver_parameters["snes_solver"]["line_search"] = "bt" # Available: basic, bt, cp, l2, nleqerr

        ### Solve the problem ###
        # self.fprint("Solving",special="header")
        # start = time.time()
        # iters, converged = nonlinear_solver.solve()
        # stop = time.time()
        # self.fprint("Total Nonlinear Iterations: {:d}".format(iters))
        # self.fprint("Converged Successfully: {0}".format(converged))
        ####################################################################

        self.fprint("Solving with {0}".format(self.nonlinear_solver))
        if self.nonlinear_solver == "newton":
            self.fprint("Relaxation parameter = {: 1.2f}".format(self.newton_relaxation))

            krylov_options = {"absolute_tolerance": 9e-3,
                              "relative_tolerance": 9e-1,
                              "maximum_iterations": 5000}

            newton_options = {"relaxation_parameter": self.newton_relaxation,
                              "maximum_iterations": 150,
                              "linear_solver": "mumps",
                              "preconditioner": "default",
                              "absolute_tolerance": 1e-6,
                              "relative_tolerance": 1e-5,
                              "error_on_nonconvergence": False,
                              "krylov_solver": krylov_options}
        
            solver_parameters = {"nonlinear_solver": "newton",
                                 "newton_solver": newton_options}

        elif self.nonlinear_solver == "snes":
            # ### Add some helper functions to solver options ###

            krylov_options = {"absolute_tolerance":  1e-12,
                              "relative_tolerance":  1e-6,
                              "maximum_iterations":  5000,
                              "monitor_convergence": True}

            solver_parameters = {"nonlinear_solver": "snes",
                                 "snes_solver": {
                                 "absolute_tolerance": 1e-6,
                                 "relative_tolerance": 1e-5,
                                 "linear_solver": "mumps",
                                 "preconditioner": "default",
                                 "maximum_iterations": 150,
                                 "error_on_nonconvergence": False,
                                 "line_search": "bt", #[basic,bt,cp,l2,nleqerr]
                                 "krylov_solver": krylov_options
                                 }}  

        else:
            raise ValueError("Unknown nonlinear solver type: {0}".format(self.nonlinear_solver))

        ### Start the Solve Process ###
        self.fprint("Solving",special="header")
        start = time.time()

        # ### Solve the Baseline Problem ###
        # solve(self.problem.F_sans_tf == 0, self.problem.up_k, self.problem.bd.bcs, solver_parameters=solver_parameters, **self.extra_kwarg)

        # ### Store the Baseline and Assign for the real solve ###
        # self.up_baseline = self.problem.up_k.copy(deepcopy=True)
        # self.problem.up_k.assign(self.up_baseline)

        ### Solve the real problem ###
        # mem0=memory_usage()[0]
        # mem_out, _ = memory_usage((solve,(self.problem.F == 0, self.problem.up_k, self.problem.bd.bcs),{"solver_parameters": solver_parameters}),max_usage=True,retval=True,max_iterations=1)
        solve(self.problem.F == 0, self.problem.up_k, self.problem.bd.bcs, solver_parameters=solver_parameters)
        stop = time.time()

        self.fprint("Solve Complete: {:1.2f} s".format(stop-start),special="footer")
        # self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))
        # self.u_k,self.p_k = self.problem.up_k.split(True)
        self.problem.u_k,self.problem.p_k = self.problem.up_k.split()

        # try:
        #     print(f"u(0, 0,150):   {self.problem.u_k([0.0, 0.0,150.0])}")
        # except:
        #     pass
        # try:
        #     print(f"u(0, 0,210):   {self.problem.u_k([0.0, 0.0,210.0])}")
        # except:
        #     pass
        # try:
        #     print(f"u(0,60,150):   {self.problem.u_k([0.0,60.0,150.0])}")
        # except:
        #     pass
        # print(f"max(u):        {self.problem.u_k.vector().max()}")
        # print(f"min(u):        {self.problem.u_k.vector().min()}")
        # print(f"integral(u_x): {assemble(self.problem.u_k[0]*dx)}")


        ### Hack into doflin adjoint to update the local controls at the start of the adjoint solve ###
        self.nu_T = project(self.problem.nu_T,self.problem.fs.Q,solver_type='gmres',preconditioner_type="hypre_amg",**self.extra_kwarg)
        if self.problem.dom.dim == 3:
            self.fprint("")
            self.fprint("Projecting Reynolds Stress")
            self.ReyStress = project(self.problem.ReyStress,self.problem.fs.T,solver_type='gmres',preconditioner_type="hypre_amg",**self.extra_kwarg)
            self.vertKE = project(self.problem.vertKE,self.problem.fs.Q,solver_type='gmres',preconditioner_type="hypre_amg",**self.extra_kwarg)

        # self.nu_T = None

        ### Save solutions ###
        if "solution" in self.params.output:
            self.fprint("Saving Solution",special="header")
            self.Save(val=self.iter_val)
            self.fprint("Finished",special="footer")

        ### calculate the power for each turbine ###
        ###################################
        ### Fix how angle is transfered ###
        ###################################
        # if self.save_power:
        #     self.SavePower(((self.iter_theta-self.problem.dom.inflow_angle)))


        ### Evaluate the objectives ###
        if self.optimizing or self.save_objective:
            self.J += self.EvaluateObjective()
            # self.J += self.objective_func(self,(self.iter_theta-self.problem.dom.inflow_angle)) 
            self.J = self.control_updater(self.J, self.problem)

        if self.save_power:
            self.EvaulatePowerFunctional()

            # print(self.outflow_markers)
            # self.J += -dot(self.problem.farm.rotor_disks,self.u_k)*dx

        # self.fprint("Speed Percent of Inflow Speed")
        # ps = []
        # for i in range(6):
        #     HH = self.problem.farm.HH[0]
        #     RD = self.problem.farm.RD[0]
        #     x_val = (i+1)*RD
        #     vel = self.problem.up_k([x_val,0,HH])
        #     vel = vel[0:3]
        #     nom = np.linalg.norm(vel)
        #     perc = nom/self.problem.bd.HH_vel
        #     ps.append(perc)
        #     self.fprint("Speed Percent at ("+repr(int(x_val))+", 0, "+repr(HH)+"): "+repr(perc))
        # print(ps)

        self.DebugOutput()

# 
# ================================================================

class IterativeSteadySolver(GenericSolver):
    """
    This solver is for solving the iterative steady state problem

    Args: 
        problem (:meth:`windse.ProblemManager.GenericProblem`): a windse problem object.
    """
    def __init__(self,problem):
        super(IterativeSteadySolver, self).__init__(problem)
        self.u_k,self.p_k = split(self.problem.up_k)

    def Solve(self):
        """
        This solves the problem setup by the problem object.
        """

        def set_solver_mode(solver_type):
            if solver_type == 'steady':
                self.problem.dt_1.assign(0)
                self.problem.dt_2.assign(1)
                self.problem.dt_3.assign(1)

                sor_vel = 0.1
                sor_pr = 0.2

            elif solver_type == 'unsteady':
                delta_t = 100.0

                self.problem.dt_1.assign(1.0/delta_t)
                self.problem.dt_2.assign(1.0/delta_t)
                self.problem.dt_3.assign(delta_t)

                sor_vel = 1.0
                sor_pr = 1.0

            else:
                raise ValueError('Solver type should be "steady" or "unsteady".')

            return sor_vel, sor_pr

        def get_relaxation_param(a_min, a_max, a_c, a_cw, x):
            
            # a_min = minimum sor value
            # a_max = maximum sor value
            # a_c = crossover location
            # a_cw = crossover width

            alpha_exp = (1/a_cw)*(a_c - x)
            alpha = 1.0/(1.0+np.exp(alpha_exp))
            alpha = alpha*(a_max-a_min) + a_min
            alpha = float(alpha)
            
            return alpha

        ### Save Files before solve ###
        self.fprint("Saving Input Data",special="header")
        if "mesh" in self.params.output:
            self.problem.dom.Save(val=self.iter_val)
        if "initial_guess" in self.params.output:
            self.problem.bd.SaveInitialGuess(val=self.iter_val)
        if "height" in self.params.output and self.problem.dom.dim == 3:
            self.problem.bd.SaveHeight(val=self.iter_val)
        if "turbine_force" in self.params.output:
            self.problem.farm.SaveActuatorDisks(val=self.iter_val)
        self.fprint("Finished",special="footer")

        self.SaveTimeSeries(0)

        # Assemble left-hand side matrices
        A1 = assemble(self.problem.F1_lhs)
        A2 = assemble(self.problem.F2_lhs)
        A3 = assemble(self.problem.F3_lhs)
        A4 = assemble(self.problem.F4_lhs)
        A5 = assemble(self.problem.F5_lhs)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in self.problem.bd.bcu]
        [bc.apply(A2) for bc in self.problem.bd.bcp]
        [bc.apply(A3) for bc in self.problem.bd.bcu]
        [bc.apply(A4) for bc in self.problem.bd.bcp]

        # Assemble right-hand side vectors
        b1 = assemble(self.problem.F1_rhs)
        b2 = assemble(self.problem.F2_rhs)
        b3 = assemble(self.problem.F3_rhs)
        b4 = assemble(self.problem.F4_rhs)
        b5 = assemble(self.problem.F5_rhs)

        # Apply bounday conditions to vectors
        [bc.apply(b1) for bc in self.problem.bd.bcu]
        [bc.apply(b2) for bc in self.problem.bd.bcp]
        [bc.apply(b3) for bc in self.problem.bd.bcu]
        [bc.apply(b4) for bc in self.problem.bd.bcp]

        # Create solvers using the set_operator method
        # which ensures preconditioners are reused when possible
        # solver_1 = PETScKrylovSolver('gmres', 'jacobi')

        vel_sol_options = ['gmres', 'none']
        pr_sol_options = ['bicgstab', 'none']

        solver_1 = PETScKrylovSolver(vel_sol_options[0], vel_sol_options[1])
        solver_1.set_operator(A1)

        solver_2 = PETScKrylovSolver(pr_sol_options[0], pr_sol_options[1])
        solver_2.set_operator(A2)

        solver_3 = PETScKrylovSolver(vel_sol_options[0], vel_sol_options[1])
        solver_3.set_operator(A3)

        solver_4 = PETScKrylovSolver(pr_sol_options[0], pr_sol_options[1])
        solver_4.set_operator(A4)

        solver_5 = PETScKrylovSolver('cg', 'jacobi')
        solver_5.set_operator(A5)



        outer_tic = time.time()


        solver_type = 'steady'
        sor_vel, sor_pr = set_solver_mode(solver_type)

        seed_velocity = False
        seed_pressure = False


        if seed_velocity:
            self.problem.u_k.assign(self.problem.bd.bc_velocity)
            self.problem.u_k_old.assign(self.problem.bd.bc_velocity)

        if seed_pressure:
            if self.params.rank == 0:
                print('Calculating initial pressure from velocity guess.')

            A1 = assemble(self.problem.F1_lhs, tensor=A1)
            [bc.apply(A1) for bc in self.problem.bd.bcu]
            solver_1.set_operator(A1)

            # Step 1: Tentative velocity step
            b1 = assemble(self.problem.F1_rhs, tensor=b1)
            [bc.apply(b1) for bc in self.problem.bd.bcu]
            solver_1.solve(self.problem.u_hat.vector(), b1)

            # Step 2: Pressure correction step
            b2 = assemble(self.problem.F2_rhs, tensor=b2)
            [bc.apply(b2) for bc in self.problem.bd.bcp]

            solver_2.solve(self.problem.p_k.vector(), b2)



        max_iter = 10

        for iter_num in range(max_iter):
            tic = time.time()

            # sor_vel = get_relaxation_param(1e-3, 0.5, 100, 10, iter_num)
            sor_vel = get_relaxation_param(1e-1, 0.5, 20, 5, iter_num)
            # sor_vel = 0.1
            sor_pr = 1.6*sor_vel

            use_simpler = False

            # if iter_num > 0.9*max_iter and solver_type == 'unsteady':
            #     solver_type = 'steady'
                # dt_1, dt_2, dt_3, sor_vel, sor_pr = set_solver_mode(solver_type)

            if use_simpler:
                A1 = assemble(self.problem.F1_lhs, tensor=A1)
                [bc.apply(A1) for bc in self.problem.bd.bcu]
                solver_1.set_operator(A1)

                # Step 1: Tentative velocity step
                b1 = assemble(self.problem.F1_rhs, tensor=b1)
                [bc.apply(b1) for bc in self.problem.bd.bcu]
                solver_1.solve(self.problem.u_hat.vector(), b1)

                # u_hat.assign((1.0-sor_vel)*u_k + sor_vel*u_hat)

                # Step 2: Pressure correction step
                b2 = assemble(self.problem.F2_rhs, tensor=b2)
                [bc.apply(b2) for bc in self.problem.bd.bcp]

                solver_2.solve(self.problem.p_k.vector(), b2)


            # Re-solve for velocity
            A3 = assemble(self.problem.F3_lhs, tensor=A3)
            [bc.apply(A3) for bc in self.problem.bd.bcu]
            solver_3.set_operator(A3)

            # Step 1: Tentative velocity step
            b3 = assemble(self.problem.F3_rhs, tensor=b3)
            [bc.apply(b3) for bc in self.problem.bd.bcu]
            solver_3.solve(self.problem.u_s.vector(), b3)

            self.problem.u_s.assign((1.0-sor_vel)*self.problem.u_k + sor_vel*self.problem.u_s)

            # Step 2: Pressure correction step
            # Don't need A4 rebuild
            b4 = assemble(self.problem.F4_rhs, tensor=b4)
            [bc.apply(b4) for bc in self.problem.bd.bcp]
            solver_4.solve(self.problem.dp.vector(), b4)


            # Step 3: Velocity correction step
            b5 = assemble(self.problem.F5_rhs, tensor=b5)
            # [bc.apply(b3) for bc in bcu]
            solver_5.solve(self.problem.du.vector(), b5)


            norm_du = norm(self.problem.du)
            norm_dp = norm(self.problem.dp)
            # u_conv.append(norm_du)
            # p_conv.append(norm_dp)

            self.problem.u_k_old.assign(self.problem.u_k)
            self.problem.p_k_old.assign(self.problem.p_k)

            self.problem.u_k.assign(self.problem.u_s + self.problem.du)

            # p_k.assign(p_k + dp)
            # p_k.assign(p_k + sor_pr*dp)

            if not use_simpler:
                # p_k.assign(p_k + dp)
                self.problem.p_k.assign(self.problem.p_k + sor_pr*self.problem.dp)

            # u_k1.assign(u_k)
            # p_k1.assign(p_k)

            if (iter_num+1) % 1 == 0:
                self.SaveTimeSeries(iter_num+1)

            toc = time.time()

            if self.params.rank == 0:
                print('Step %3d of %3d, sor = (%.3f, %.3f): | du | = %.6e, | dp | = %.6e, CPU Time = %.2f s' % (iter_num+1, max_iter, sor_vel, sor_pr, norm_du, norm_dp, toc-tic))

        outer_toc = time.time()
        if self.params.rank == 0:
            print('TOTAL CPU TIME = %f seconds' % (outer_toc - outer_tic))

        ### Evaluate the objectives ###
        if self.optimizing or self.save_objective:
            self.J += self.EvaluateObjective()
            # self.J += self.objective_func(self,(self.iter_theta-self.problem.dom.inflow_angle)) 
            self.J = self.control_updater(self.J, self.problem)

        if self.save_power and "power":
            self.EvaulatePowerFunctional()

    @no_annotations
    def SaveTimeSeries(self, simTime):
        # if hasattr(self.problem,"tf_save"):
        #     self.problem.tf_save.vector()[:] = 0
        #     for fun in self.problem.tf_list:
        #         self.problem.tf_save.vector()[:] = self.problem.tf_save.vector()[:] + fun.vector()[:]
        # else:
        #     self.problem.tf_save = Function(self.problem.fs.V)
        #     for fun in self.problem.tf_list:
        #         self.problem.tf_save.vector()[:] = self.problem.tf_save.vector()[:] + fun.vector()[:]


        if self.first_save:
            self.velocity_file = self.params.Save(self.problem.u_k,"velocity",subfolder="timeSeries/",val=simTime)
            self.pressure_file   = self.params.Save(self.problem.p_k,"pressure",subfolder="timeSeries/",val=simTime)
            # self.turb_force_file   = self.params.Save(self.problem.tf_save,"turbine_force",subfolder="timeSeries/",val=simTime)
            # if self.optimizing:
                # self.adj_tape_file = XDMFFile(self.params.folder+"timeSeries/global_adjoint.xdmf")
                # self.problem.u_k1.assign(Marker(self.problem.u_k,simTime,self.adj_tape_file))
            self.first_save = False
        else:
            self.params.Save(self.problem.u_k,"velocity",subfolder="timeSeries/",val=simTime,file=self.velocity_file)
            self.params.Save(self.problem.p_k,"pressure",subfolder="timeSeries/",val=simTime,file=self.pressure_file)
            # self.params.Save(self.problem.tf_save,"turbine_force",subfolder="timeSeries/",val=simTime,file=self.turb_force_file)

# ================================================================

class UnsteadySolver(GenericSolver):
    """
    This solver is for solving an unsteady problem.  As such, it contains
    additional time-stepping features and functions not present in other solvers.
    This solver can only be used if an unsteady problem has been specified in
    the input file.

    Args: 
        problem (:meth:`windse.ProblemManager.GenericProblem`): a windse problem object.
    """
    def __init__(self,problem):
        super(UnsteadySolver, self).__init__(problem)

    # ================================================================

    def Solve(self):
        # Start the unsteady solver ONLY if an unsteady problem has been created
        if self.problem.params["problem"]["type"] == 'unsteady':
            self.fprint("Solving with UnsteadySolver", special="header")
        else:
            raise ValueError("UnsteadySolver can only be run with ProblemType = unsteady, not %s" \
                % (self.problem.params["problem"]["type"]))

        if type(self.final_time) == str and self.final_time.lower() == 'none':
            self.pseudo_steady = True
            self.fprint('Found option "None" for final_time, ')
            self.fprint('Running until unsteady solver is converged.')
            self.final_time = 1000000.0
            #self.final_time = 10.0

        # ================================================================
        
        # Start a counter for the total simulation time
        self.fprint("dt: %.4f" % (self.problem.dt))
        self.fprint("Final Time: %.1f" % (self.final_time))

        # Calculate the time to start recording if optimizing
        if self.optimizing:
            if self.record_time == "computed":
                self.record_time = 1.0*(self.wake_RD*self.problem.farm.RD[0]/(self.problem.bd.HH_vel*0.75))
                if self.final_time < self.record_time + 60.0/self.problem.rpm:
                    self.fprint("Warning: Final time is too small... overriding")
                    self.final_time = self.record_time + 60.0/self.problem.rpm
                    self.fprint("         New Final Time: {:1.2f} s".format(self.final_time))
            elif self.record_time == "last":
                self.record_time = self.final_time
            self.fprint("Recording Time: %.1f" % (self.record_time))
            # self.record_time = 250.0
            # self.record_time = self.final_time-2.0*saveInterval
        else:
            self.record_time = 0.0
        self.problem.record_time = self.record_time

        # Calculate what time to start the averaging process
        init_average = True
        if self.pseudo_steady:
            two_flowthrough_time = 2.0*(self.problem.dom.x_range[1] - self.problem.dom.x_range[0])/self.params['boundary_conditions']['HH_vel']
            average_start_time = two_flowthrough_time
            self.fprint('Start averaging after two flow-throughs, or %.1f seconds' % (two_flowthrough_time))
        else:
            average_start_time = 5.0

        # ================================================================

        # Start a counter for the number of saved files
        saveCount = 0
        save_next_timestep = False

        # Generate file pointers for saved output
        # FIXME: This should use the .save method
        # fp = []
        # fp.append(File("%s/timeSeries/velocity.pvd" % (self.problem.dom.params.folder)))
        # fp.append(File("%s/timeSeries/pressure.pvd" % (self.problem.dom.params.folder)))
        # fp.append(File("%s/timeSeries/nu_T.pvd" % (self.problem.dom.params.folder)))

        # if "turbine_force" in self.params.output:
        #     fp.append(File("%s/timeSeries/turbineForce.pvd" % (self.problem.dom.params.folder)))

        # Save first timestep (create file pointers for first call)
        self.SaveTimeSeries(self.simTime, 0.0)

        self.fprint("Saving Input Data",special="header")
        if "mesh" in self.params.output:
            self.problem.dom.Save(val=self.iter_val)
        if "initial_guess" in self.params.output:
            self.problem.bd.SaveInitialGuess(val=self.iter_val)
        if "height" in self.params.output and self.problem.dom.dim == 3:
            self.problem.bd.SaveHeight(val=self.iter_val)
        if "turbine_force" in self.params.output:
            self.problem.farm.save_functions(val=self.iter_val)

        self.fprint("Finished",special="footer")

        # ================================================================

        self.fprint("")
        self.fprint("Calculating Boundary Conditions")

        # FIXME: This should use the boundary information in self.problem.bd.bcs
        # bcu, bcp = self.GetBoundaryConditions(0.0)
        self.problem.dom.RecomputeBoundaryMarkers(0.0)

        # ================================================================

        self.fprint("Assembling time-independent matrices")

        # Assemble left-hand side matrices
        A1 = assemble(self.problem.a1)
        A2 = assemble(self.problem.a2)
        A3 = assemble(self.problem.a3)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in self.problem.bd.bcu]
        [bc.apply(A2) for bc in self.problem.bd.bcp]

        # Assemble right-hand side vector
        b1 = assemble(self.problem.L1)
        b2 = assemble(self.problem.L2)
        b3 = assemble(self.problem.L3)

        # Apply bounday conditions to vectors
        [bc.apply(b1) for bc in self.problem.bd.bcu]
        [bc.apply(b2) for bc in self.problem.bd.bcp]

        timing = np.zeros(5)

        # sol = ['cg', 'bicgstab', 'gmres']
        # pre = ['amg', 'hypre_amg', 'hypre_euclid', 'hypre_parasails', 'jacobi', 'petsc_amg', 'sor']

        # cl_iterator = int(self.params['solver']['cl_iterator'])

        # if cl_iterator > 0:
        #     ind = cl_iterator-1
        #     sol_choice = sol[int(np.floor(ind/7))]
        #     pre_choice = pre[int(ind%7)]

        # else:
        #     sol_choice = 'gmres'
        #     pre_choice = 'default'

        solver_1 = PETScKrylovSolver('gmres', 'jacobi')
        # solver_1 = PETScKrylovSolver('gmres', 'default')
        # solver_1 = PETScKrylovSolver('default', 'default')
        solver_1.set_operator(A1)

        solver_2 = PETScKrylovSolver('gmres', 'petsc_amg')
        # solver_2 = PETScKrylovSolver('gmres', 'hypre_amg')
        # solver_2 = PETScKrylovSolver('default', 'default')
        solver_2.set_operator(A2)

        solver_3 = PETScKrylovSolver('cg', 'jacobi')
        # solver_3 = PETScKrylovSolver('gmres', 'default')
        # solver_3 = PETScKrylovSolver('default', 'default')
        solver_3.set_operator(A3)

        pr = cProfile.Profile()

        # ================================================================

        self.fprint('Turbine Parameters', special='header')
        self.fprint('Hub Height: %.1f' % (self.problem.farm.turbines[0].HH))
        self.fprint('Yaw: %.4f' % (self.problem.farm.turbines[0].yaw))
        self.fprint('Radius: %.1f' % (self.problem.farm.turbines[0].radius))
        self.fprint('', special='footer')

        self.fprint("Solving",special="header")
        self.fprint("Sim Time | Next dt | U_max")
        self.fprint("--------------------------")

        def save_adj(adj):
            # print(len(self.adj_time_list),self.adj_time_iter)
            if self.adj_time_iter < len(self.adj_time_list):
                t = np.flip(self.adj_time_list)[self.adj_time_iter]
                # print(t,t % self.save_interval)
                if t % self.save_interval == 0:
                    adj.rename("adjoint","adjoint")
                    self.adj_file.write(adj,t)
            self.adj_time_iter += 1

        # Initialize need loop objects
        start = time.time()
        # self.problem.dt_sum = 1.0
        self.problem.dt_sum = 0.0
        J_old = 0
        J_diff_old = 100000
        min_count = 0
        simIter = 0
        stable = False

        if self.problem.farm.turbines[0].type == "line" or self.problem.farm.turbines[0].type == "dolfin_line":
            tip_speed = self.problem.farm.turbines[0].rpm*2.0*np.pi*self.problem.farm.turbines[0].radius/60.0
        else:
            tip_speed = 0.0

        # self.problem.alm_power_sum = 0.0

        while not stable and self.simTime < self.final_time:

            # add a fake block that allows us to update the control while dolfin_adjoint is doing it's thing
            # Need to allow processes which don't own the above points to fail gracefully
            # try:
            #     print(f"u(0, 0,150):   {self.problem.u_k([0.0, 0.0,150.0])}")
            # except:
            #     pass
            # try:
            #     print(f"u(0, 0,210):   {self.problem.u_k([0.0, 0.0,210.0])}")
            # except:
            #     pass
            # try:
            #     print(f"u(0,60,150):   {self.problem.u_k([0.0,60.0,150.0])}")
            # except:
            #     pass
            # print(f"max(u):        {self.problem.u_k.vector().max()}")
            # print(f"min(u):        {self.problem.u_k.vector().min()}")
            # print(f"integral(u_x): {assemble(self.problem.u_k[0]*dx)}")
            self.J = self.control_updater(self.J, self.problem, time=self.simTime)

            self.problem.bd.UpdateVelocity(self.simTime)

            # Record the "old" max velocity (before this update)
            u_max_k1 = max(tip_speed, self.problem.u_k.vector().max())

            # Step 1: Tentative velocity step
            tic = time.time()
            # solve(self.problem.a1==self.problem.L1, self.problem.u_k, bcs=self.problem.bd.bcu)
            # solve(self.problem.a1==self.problem.L1, self.problem.u_k, bcs=self.problem.bd.bcu, solver_parameters={"linear_solver": "gmres","preconditioner": "jacobi"})

            b1 = assemble(self.problem.L1, tensor=b1)
            [bc.apply(b1) for bc in self.problem.bd.bcu]
            if self.optimizing:
                # solve(A1, self.problem.u_k.vector(), b1, 'gmres', 'default',adj_cb=save_adj)
                # solve(A1, self.problem.u_k.vector(), b1, 'gmres', 'default')
                solver_1.solve(self.problem.u_k.vector(), b1)
            else:
                # solve(A1, self.problem.u_k.vector(), b1, 'gmres', 'default')
                solver_1.solve(self.problem.u_k.vector(), b1)
            # print("uk("+repr(self.simTime)+")   = "+repr(np.mean(self.problem.u_k.vector()[:])))
            # print("assemble(func*dx): " + repr(float(assemble(inner(self.problem.u_k,self.problem.u_k)*dx))))
            toc = time.time()
            timing[0] += toc - tic

            # Step 2: Pressure correction step
            tic = time.time()
            # solve(self.problem.a2==self.problem.L2, self.problem.p_k, bcs=self.problem.bd.bcp)
            # solve(self.problem.a2==self.problem.L2, self.problem.p_k, bcs=self.problem.bd.bcp, solver_parameters={"linear_solver": "gmres","preconditioner": "petsc_amg"})
            
            b2 = assemble(self.problem.L2, tensor=b2)
            [bc.apply(b2) for bc in self.problem.bd.bcp]
            if self.optimizing:
                # solve(A2, self.problem.p_k.vector(), b2, 'gmres', 'hypre_amg')
                solver_2.solve(self.problem.p_k.vector(), b2)
            else:
                solver_2.solve(self.problem.p_k.vector(), b2)
            # print("uk("+repr(self.simTime)+")   = "+repr(np.mean(self.problem.p_k.vector()[:])))
            # print("assemble(func*dx): " + repr(float(assemble(inner(self.problem.p_k,self.problem.p_k)*dx))))
            toc = time.time()
            timing[1] += toc - tic

            # Step 3: Velocity correction step
            tic = time.time()
            # solve(self.problem.a3==self.problem.L3, self.problem.u_k)
            # solve(self.problem.a3==self.problem.L3, self.problem.u_k, solver_parameters={"linear_solver": "cg","preconditioner": "jacobi"})
            
            b3 = assemble(self.problem.L3, tensor=b3)
            if self.optimizing:
                # solve(A3, self.problem.u_k.vector(), b3, 'gmres', 'default')
                solver_3.solve(self.problem.u_k.vector(), b3)
            else:
                solver_3.solve(self.problem.u_k.vector(), b3)
            # print("uk("+repr(self.simTime)+")   = "+repr(np.mean(self.problem.u_k.vector()[:])))
            # print("assemble(func*dx): " + repr(float(assemble(inner(self.problem.u_k,self.problem.u_k)*dx))))
            toc = time.time()
            timing[2] += toc - tic

            # Old <- New update step
            self.problem.u_k2.assign(self.problem.u_k1)
            self.problem.u_k1.assign(self.problem.u_k)
            self.problem.p_k1.assign(self.problem.p_k)

            # Record the updated max velocity
            u_max = max(tip_speed, self.problem.u_k.vector().max())

            # Update the simulation time
            self.simTime_prev = self.simTime
            self.simTime += self.problem.dt

            # Compute Reynolds Stress
            if 'KE_entrainment' in self.objective_type.keys():
                if self.simTime >= self.u_avg_time and self.simTime < self.record_time:
                    self.problem.uk_sum.assign(self.problem.uk_sum+self.problem.dt_c*self.problem.u_k)
                elif self.simTime >= self.record_time:
                    self.problem.vertKE = (self.problem.u_k[0]-self.problem.uk_sum[0]/(self.record_time-self.u_avg_time))*(self.problem.u_k[2]-self.problem.uk_sum[2]/(self.record_time-self.u_avg_time))*(self.problem.uk_sum[0]/(self.record_time-self.u_avg_time))
            
            if self.save_all_timesteps:
                self.SaveTimeSeries(self.simTime,simIter)
            elif save_next_timestep:
                # Read in new inlet values
                # bcu = self.updateInletVelocityFromFile(saveCount, bcu)
                
                # Clean up self.simTime to avoid accumulating round-off error
                saveCount += 1
                self.simTime = self.save_interval*saveCount

                # Save output files
                # self.SaveTimeSeries(fp, self.simTime)
                self.SaveTimeSeries(self.simTime,simIter)

            # Adjust the timestep size, dt, for a balance of simulation speed and stability
            save_next_timestep = self.AdjustTimestepSize(save_next_timestep, self.save_interval, self.simTime, u_max, u_max_k1)

            # Update the turbine force
            tic = time.time()
            if self.problem.farm.turbines[0].type == "line" or self.problem.farm.turbines[0].type == "dolfin_line":

                self.problem.alm_power = 0.0
                # pass

                # t1 = time.time()
                pr.enable()
                if 'dolfin' in self.problem.farm.turbines[0].type:
                    self.problem.ComputeTurbineForceTerm(self.problem.u_k, TestFunction(self.problem.fs.V), self.problem.dom.inflow_angle, simTime=self.simTime, simTime_prev=self.simTime_prev, dt=self.problem.dt)
                else:
                    # updated method from dev, is this necessary? this may not work if there are self.farm.too_many_turbines
                    self.problem.farm.update_turbine_force(self.problem.u_k, self.problem.dom.inflow_angle, self.problem.fs, simTime=self.simTime, simTime_prev=self.simTime_prev, dt=self.problem.dt)
                    # new_tf = self.problem.ComputeTurbineForce(self.problem.u_k, self.problem.dom.inflow_angle, simTime=self.simTime, simTime_prev=self.simTime_prev, dt=self.problem.dt)
                    # self.problem.tf.assign(new_tf)
                pr.disable()

                # # t2 = time.time()
                # # print(t2-t1)

                # Power [=] N*m*rads/s 
                # self.problem.alm_power = self.problem.rotor_torque*(2.0*np.pi*self.problem.rpm/60.0)
                # self.problem.alm_power_dolfin = self.problem.rotor_torque_dolfin*(2.0*np.pi*self.problem.rpm/60.0)
                
                # # self.problem.alm_power_sum += self.problem.alm_power*self.problem.dt
                # # self.problem.alm_power_average = self.problem.alm_power_sum/self.simTime

                # self.problem.alm_power_dolfin = self.problem.farm.compute_power(self.problem.u_k, self.problem.dom.inflow_angle)
                # output_str = 'Rotor Power (dolfin, solver): %s MW' % (self.problem.alm_power_dolfin/1.0e6)
                # self.fprint(output_str)

            else:
                # This is a hack to avoid errors when using something other than ALM
                # e.g., actuator disks
                self.problem.alm_power = np.zeros(self.problem.farm.numturbs)
                self.problem.alm_power_dolfin = np.zeros(self.problem.farm.numturbs)

                # exit()
            toc = time.time()
            timing[3] += toc - tic

            if self.simTime > average_start_time:
                if init_average:

                    convergence_history = []

                    average_vel_sum = Function(self.problem.fs.V)
                    average_vel_1 = Function(self.problem.fs.V)
                    average_vel_2 = Function(self.problem.fs.V)

                    average_vel_sum.vector()[:] = self.problem.u_k1.vector()[:]*self.problem.dt_previous
                    average_vel_1.vector()[:] = average_vel_sum.vector()[:]/(self.simTime-average_start_time)
                    average_vel_2.vector()[:] = 0.0

                    average_power = self.problem.alm_power*self.problem.dt

                    init_average = False
                else:
                    average_vel_sum.vector()[:] += self.problem.u_k1.vector()[:]*self.problem.dt_previous

                    average_vel_2.vector()[:] = average_vel_1.vector()[:]
                    average_vel_1.vector()[:] = average_vel_sum.vector()[:]/(self.simTime-average_start_time)

                    norm_diff = norm(average_vel_2.vector() - average_vel_1.vector(), 'linf')
                    self.fprint('Change Between Steps (norm)'+ repr(norm_diff))


                    if norm_diff < 1e-4:
                        # If the average is no longer changing, stop this run
                        convergence_history.append(1)

                        # Only mark "stable" if the last 2 tests have passed
                        if convergence_history[-1] == 1 and convergence_history[-2] == 1:
                            stable = True
                    else:
                        convergence_history.append(0)

                    average_power += self.problem.alm_power*self.problem.dt


            # # Adjust the timestep size, dt, for a balance of simulation speed and stability
            # save_next_timestep = self.AdjustTimestepSize(save_next_timestep, self.save_interval, self.simTime, u_max, u_max_k1)

            if self.save_objective or (self.optimizing and self.simTime >= self.record_time and not self.pseudo_steady):
                J_next = self.EvaluateObjective()

            # Calculate the objective function
            if self.optimizing and self.simTime >= self.record_time and not self.pseudo_steady:

                # Append the current time step for post production
                self.adj_time_list.append(self.simTime)

                self.J += float(self.problem.dt)*J_next
                self.problem.dt_sum += self.problem.dt 
                J_new = float(self.J/self.problem.dt_sum)

                ### TODO, replace this with an actual stabilization criteria such as relative difference tolerance
                # Check the change in J with respect to time and check if we are "stable" i.e. hit the required number of minimums
                J_diff = J_new-J_old
                # if J_diff_old <= 0 and J_diff > 0 and self.min_total:
                #     if min_count == self.min_total:
                #         stable = True
                #     else:
                #         min_count += 1
                # J_diff_old = J_diff
                J_old = J_new

                # Another stable checking method that just looks at the difference
                # if abs(J_diff) <= 0.001:
                #     stable = True

                self.fprint("Current Objective Value: "+repr(float(self.J/self.problem.dt_sum)))
                self.fprint("Change in Objective    : "+repr(float(J_diff)))



            # to only call the power functional once, check if a) the objective is the power, b) that we are before record time
            if self.save_power:
                self.EvaulatePowerFunctional()

            # After changing timestep size, A1 must be reassembled
            # FIXME: This may be unnecessary (or could be sped up by changing only the minimum amount necessary)

            tic = time.time()
            A1 = assemble(self.problem.a1, tensor=A1)
            [bc.apply(A1) for bc in self.problem.bd.bcu]
            solver_1.set_operator(A1)
            toc = time.time()
            timing[4] += toc - tic

            # Print some solver statistics
            self.fprint("%8.2f | %7.2f | %5.2f" % (self.simTime, self.problem.dt, u_max))
            simIter+=1

        if self.pseudo_steady:
            self.J = self.EvaluateObjective()

        elif (self.optimizing or self.save_objective):
            # if self.problem.dt_sum > 0.0:
            self.J = self.J/float(self.problem.dt_sum)


        if self.simTime > average_start_time:
            average_vel_sum.vector()[:] = average_vel_sum.vector()[:]/(self.simTime-average_start_time)
            # fp = File('./output/%s/average_vel_sum.pvd' % (self.params.name))
            fp = File(f"{self.params.folder}timeSeries/average_vel_sum.pvd")
            average_vel_sum.rename('average_vel_sum', 'average_vel_sum')
            fp << average_vel_sum

            average_power = average_power/(self.simTime-average_start_time)
            # self.fprint('AVERAGE Rotor Power: %.6f MW' % (average_power/1e6))
            try:
                output_str = 'AVERAGE Rotor Power: %s MW' % (np.array2string(average_power/1.0e6, precision=9, separator=', '))
            except:
                output_str = 'AVERAGE Rotor Power: %s MW' % (average_power)
            self.fprint(output_str)

        # add a fake block that allows us to update the control while dolfin_adjoint is doing it's thing
        # try:
        #     print(f"u(0, 0,150):   {self.problem.u_k([0.0, 0.0,150.0])}")
        # except:
        #     pass
        # try:
        #     print(f"u(0, 0,210):   {self.problem.u_k([0.0, 0.0,210.0])}")
        # except:
        #     pass
        # try:
        #     print(f"u(0,60,150):   {self.problem.u_k([0.0,60.0,150.0])}")
        # except:
        #     pass

        # print(f"max(u):        {self.problem.u_k.vector().max()}")
        # print(f"min(u):        {self.problem.u_k.vector().min()}")
        # print(f"integral(u_x): {assemble(self.problem.u_k[0]*dx)}")
        self.J = self.control_updater(self.J, self.problem, time=self.simTime)
        stop = time.time()

        self.fprint('================================================================')
        # self.fprint('Solver:              %s' % (sol_choice))
        # self.fprint('Preconditioner:      %s' % (pre_choice))
        total_time = stop - start
        self.fprint('Assembling A1:       %.2f sec (%.1f%%)' % (timing[4], 100.0*timing[4]/total_time))
        self.fprint('Tentative Velocity:  %.2f sec (%.1f%%)' % (timing[0], 100.0*timing[0]/total_time))
        self.fprint('Pressure Correction: %.2f sec (%.1f%%)' % (timing[1], 100.0*timing[1]/total_time))
        self.fprint('Velocity Update:     %.2f sec (%.1f%%)' % (timing[2], 100.0*timing[2]/total_time))
        self.fprint('ALM Calculation:     %.2f sec (%.1f%%)' % (timing[3], 100.0*timing[3]/total_time))
        self.fprint('================================================================')
        self.fprint("Finished",special="footer")
        self.fprint("Solve Complete: {:1.2f} s".format(stop-start),special="footer")

        # if self.params.rank == 0:
        #     # pr.print_stats(25, sort = 'time')
        #     ps = pstats.Stats(pr).strip_dirs().sort_stats('cumulative')
        #     ps.print_stats(50)
        #     pr.dump_stats('%s/profiling_alm.prof' % (self.params.folder))


    # ================================================================
    @no_annotations
    def SaveTimeSeries(self, simTime, simIter=None):

        self.DebugOutput(simTime,simIter)
        ### TODO THIS NEED TO BE CLEAN TO ACCOUNT FOR DISKS

        if self.problem.farm.turbines[0].type == "line":
            if hasattr(self.problem,"tf_save"):
                self.problem.tf_save.vector()[:] = 0
                for fun in self.problem.farm.turbines:
                    self.problem.tf_save.vector()[:] = self.problem.tf_save.vector()[:] + fun.tf.vector()[:]
            else:
                self.problem.tf_save = Function(self.problem.fs.V)
                for fun in self.problem.farm.turbines:
                    self.problem.tf_save.vector()[:] = self.problem.tf_save.vector()[:] + fun.tf.vector()[:]
        else:
            if hasattr(self.problem,"tf_save"):
                self.problem.tf_save = project(sum(self.problem.farm.tf_list), self.problem.fs.V, solver_type='cg')
            else:
                self.problem.tf_save = Function(self.problem.fs.V)
                self.problem.tf_save = project(sum(self.problem.farm.tf_list), self.problem.fs.V, solver_type='cg')


        if self.first_save:
            self.velocity_file = self.params.Save(self.problem.u_k,"velocity",subfolder="timeSeries/",val=simTime)
            self.pressure_file   = self.params.Save(self.problem.p_k,"pressure",subfolder="timeSeries/",val=simTime)
            self.turb_force_file   = self.params.Save(self.problem.tf_save,"turbine_force",subfolder="timeSeries/",val=simTime)
            # if self.optimizing:
                # self.adj_tape_file = XDMFFile(self.params.folder+"timeSeries/global_adjoint.xdmf")
                # self.problem.u_k1.assign(Marker(self.problem.u_k,simTime,self.adj_tape_file))
            self.first_save = False
        else:
            self.params.Save(self.problem.u_k,"velocity",subfolder="timeSeries/",val=simTime,file=self.velocity_file)
            self.params.Save(self.problem.p_k,"pressure",subfolder="timeSeries/",val=simTime,file=self.pressure_file)
            self.params.Save(self.problem.tf_save,"turbine_force",subfolder="timeSeries/",val=simTime,file=self.turb_force_file)
            # if self.optimizing:
                # self.problem.u_k1.assign(Marker(self.problem.u_k,simTime,self.adj_tape_file))

        # # Save velocity files (pointer in fp[0])
        # self.problem.u_k.rename('Velocity', 'Velocity')
        # fp[0] << (self.problem.u_k, simTime)

        # # Save pressure files (pointer in fp[1])
        # self.problem.p_k.rename('Pressure', 'Pressure')
        # fp[1] << (self.problem.p_k, simTime)

        # # Save eddy viscosity files (pointer in fp[2])
        # # nu_T_val = project(self.problem.nu_T, self.problem.fs.Q, solver_type='gmres')
        # # nu_T_val.rename('nu_T', 'nu_T')
        # # fp[2] << (nu_T_val, simTime)

        # # Save turbine force files (pointer in fp[3])
        # # if "turbine_force" in self.params.output:

        # workaround = False

        # if workaround:
        #     tf_value = project(self.problem.tf, self.problem.fs.V, solver_type='gmres')
        #     tf_value.rename('Turbine_Force', 'Turbine_Force')
        #     fp[3] << (tf_value, simTime)
        # else:
        #     self.problem.tf.rename('Turbine_Force', 'Turbine_Force')
        #     fp[3] << (self.problem.tf, simTime)


    # ================================================================

    def AdjustTimestepSize(self, save_next_timestep, saveInterval, simTime, u_max, u_max_k1):
        # Save the old dt for reference
        self.problem.dt_previous = self.problem.dt

        # Set the CFL target (0.2 is a good value for stability and speed, YMMV)
        # cfl_target = 0.5
        cfl_target = 1.0
        cfl_target = float(self.problem.params["solver"]["cfl_target"])

        # Enforce a minimum timestep size
        dt_min = 0.01

        # Calculate the change in velocity using a first-order, backward difference
        dudt = u_max - u_max_k1

        # Calculate the projected velocity
        u_max_projected = u_max + dudt

        # Calculate the ideal timestep size (ignore file output considerations for now)
        dt_new = cfl_target * self.problem.dom.global_hmin / u_max_projected

        # Move to larger dt slowly (smaller dt happens instantly)
        if dt_new > self.problem.dt:
            # Amount of new dt to use: 0 = none, 1 = all
            SOR = 0.5
            dt_new = SOR*dt_new + (1.0-SOR)*self.problem.dt

        # dt_new = dt_new/2.0

        # Calculate the time remaining until the next file output
        time_remaining = saveInterval - (simTime % saveInterval)

        # If the new timestep would jump past a save point, modify the new timestep size
        if not save_next_timestep and dt_new + dt_min >= time_remaining:
            dt_new = time_remaining
            save_next_timestep = True
        else:
            save_next_timestep = False

        # dt_new = 0.04
        if self.params.num_procs > 1:
            dt_new_vec = np.zeros(self.params.num_procs, dtype=np.float64)
            dt_new = np.array(dt_new, dtype=np.float64)
            self.params.comm.Allgather(dt_new, dt_new_vec)
            dt_new = np.amin(dt_new_vec)

        # print('Rank %d dt = %.15e' % (self.params.rank, dt_new))

        # Update both the Python variable and FEniCS constant
        self.problem.dt = dt_new
        self.problem.dt_c.assign(dt_new)


        # float(self.problem.dt_c) # to get the regular ol' variable

        return save_next_timestep

    # ================================================================

    def UpdateActuatorLineForce(self, simTime):

        def rot_x(theta):
            Rx = np.array([[1, 0, 0],
                           [0, np.cos(theta), -np.sin(theta)],
                           [0, np.sin(theta), np.cos(theta)]])

            return Rx

        def rot_y(theta):
            Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                           [0, 1, 0],
                           [-np.sin(theta), 0, np.cos(theta)]])
            
            return Ry

        def rot_z(theta):
            Rz = np.array([[np.cos(theta), -np.sin(theta), 0],
                           [np.sin(theta), np.cos(theta), 0],
                           [0, 0, 1]])
            
            return Rz

        #================================================================
        # Get Mesh Properties
        #================================================================

        ndim = self.problem.dom.dim

        # Get the coordinates of the vector function space
        coords = self.problem.fs.V.tabulate_dof_coordinates()
        coords = np.copy(coords[0::self.problem.dom.dim, :])


        # Resape a linear copy of the coordinates for every mesh point
        coordsLinear = np.copy(coords.reshape(-1, 1))

        #================================================================
        # Set Turbine and Fluid Properties
        #================================================================

        # Set the density
        rho = 1.0

        # Set the hub height
        hub_height = self.problem.farm.HH[0] # For a SWIFT turbine

        # Get the hub-height velocity
        u_inf = 8.0

        # Set the rotational speed of the turbine
        RPM = 15.0

        # Set the yaw of the turbine
        yaw = self.problem.farm.yaw[0]

        # Set the number of blades in the turbine
        num_blades = 3

        # Blade length (turbine radius)
        L = self.problem.farm.radius[0] # For a SWIFT turbine 27 m in diameter

        # Chord length
        c = L/20.0

        # Width of Gaussian
        # eps = 2.5*c
        eps = 2.0*self.problem.dom.mesh.hmin()/np.sqrt(3)

        # print(self.problem.farm.x)
        # print(self.problem.farm.y)
        # print(self.problem.farm.HH)
        # print(self.problem.farm.yaw)
        # print(self.problem.farm.RD)
        # print(self.problem.farm.radius)

        # Discretize each blade into separate nodes
        num_actuator_nodes = 10

        #================================================================
        # Set Derived Constants
        #================================================================

        # Calculate the blade velocity
        period = 60.0/RPM
        tip_speed = np.pi*2.0*L*RPM/60.0
        blade_vel = np.vstack((np.zeros(num_actuator_nodes),
                               np.zeros(num_actuator_nodes),
                               np.linspace(0.0, tip_speed, num_actuator_nodes)))

        # Set the initial angle of each blade
        theta_vec = np.linspace(0.0, 2.0*np.pi, num_blades+1)
        theta_vec = theta_vec[0:num_blades]

        # Calculate discrete node positions
        rdim = np.linspace(0.0, L, num_actuator_nodes)

        # Calculate width of individual blade segment
        w = rdim[1] - rdim[0]

        # Calculate an array describing the x, y, z position of each point
        xblade = np.vstack((np.zeros(num_actuator_nodes),
                            rdim,
                            np.zeros(num_actuator_nodes)))

        #================================================================
        # Begin Calculating Turbine Forces
        #================================================================

        # Lift and drag coefficient (could be an array and you interpolate the value based on R)
        # cl_dolf = Constant((np.linspace(1.5, 0.5, num_actuator_nodes)))
        # cd_dolf = Constant((np.ones(num_actuator_nodes)))
        # cl = cl_dolf.values()
        # cd = cd_dolf.values()

        cl = np.linspace(0.0, 2.0, num_actuator_nodes) # Uncomment for controllability study
        cd = np.linspace(2.0, 0.0, num_actuator_nodes)
        # cl = np.linspace(2.0, 0.0, num_actuator_nodes) # Uncomment for controllability study
        # cd = np.linspace(0.0, 2.0, num_actuator_nodes)
        # cl = np.ones(num_actuator_nodes)
        # cd = np.ones(num_actuator_nodes)


        # Create space to hold the vector values
        tf_vec = np.zeros(np.size(coords))
        lift_force = np.zeros((np.shape(coords)[0], ndim))
        drag_force = np.zeros((np.shape(coords)[0], ndim))

        # tf_lift_vec = np.zeros(np.size(coords))
        # tf_drag_vec = np.zeros(np.size(coords))

        # Calculate the blade position based on current simTime and turbine RPM
        theta_offset = simTime/period*2.0*np.pi

        # Treat each blade separately
        for theta_0 in theta_vec:
            theta = theta_0 + theta_offset

            # Generate a rotation matrix for this turbine blade
            Rx = rot_x(theta)
            Rz = rot_z(yaw)

            # Rotate the entire [x; y; z] matrix using this matrix, then shift to the hub height
            xblade_rotated = np.dot(Rz, np.dot(Rx, xblade))
            xblade_rotated[2, :] += hub_height

            # Tile the blade coordinates for every mesh point, [numGridPts*ndim x num_actuator_nodes]
            xblade_rotated_full = np.tile(xblade_rotated, (np.shape(coords)[0], 1))

            # Subtract and square to get the dx^2 values in the x, y, and z directions
            dx_full = (coordsLinear - xblade_rotated_full)**2

            # Add together to get |x^2 + y^2 + z^2|^2
            dist2 = dx_full[0::ndim] + dx_full[1::ndim] + dx_full[2::ndim]

            # Set if using local velocity around inidividual nodes
            using_local_velocity = False
        
            if using_local_velocity:
                # Generate the fluid velocity from the actual node locations in the flow
                u_fluid = np.zeros((3, num_actuator_nodes))
                
                for k in range(num_actuator_nodes):
                    u_fluid[:, k] = self.problem.u_k1(xblade_rotated[0, k],
                                                      xblade_rotated[1, k],
                                                      xblade_rotated[2, k])
                                    
            else:
                # Generate the fluid velocity analytically using the hub height velocity
                # u_inf_vec = u_inf*np.ones(num_actuator_nodes)
                
                # u_fluid = np.vstack((u_inf_vec,
                #                      np.zeros(num_actuator_nodes),
                #                      np.zeros(num_actuator_nodes)))
                u_fluid = np.zeros((3, num_actuator_nodes))
                
                for k in range(num_actuator_nodes):
                    u_fluid[0, k] = 8.0*(xblade_rotated[2, k]/hub_height)**0.18

            
            # Rotate the blade velocity in the global x, y, z, coordinate system
            blade_vel_rotated = np.dot(Rz, np.dot(Rx, -blade_vel))
                            
            # Form the total relative velocity vector (including velocity from rotating blade)
            u_rel = u_fluid + blade_vel_rotated
            
            # Create unit vectors in the direction of u_rel
            u_rel_mag = np.linalg.norm(u_rel, axis=0)
            u_rel_mag[u_rel_mag < 1e-6] = 1e-6
            u_unit = u_rel/u_rel_mag
            
            # Calculate the lift and drag forces using the relative velocity magnitude
            lift = (0.5*cl*rho*c*w*u_rel_mag**2)/(eps**3 * np.pi**1.5)
            drag = (0.5*cd*rho*c*w*u_rel_mag**2)/(eps**3 * np.pi**1.5)
            
            # Calculate the force at every mesh point due to every node [numGridPts x NumActuators]
            nodal_lift = lift*np.exp(-dist2/eps**2)
            nodal_drag = drag*np.exp(-dist2/eps**2)
            
            # Calculate a vector in the direction of the blade
            blade_unit = xblade_rotated[:, -1] - np.array([0.0, 0.0, hub_height])  
            
            for k in range(num_actuator_nodes):
                # The drag unit simply points opposite the relative velocity unit vector
                drag_unit = -u_unit[:, k]
                
                # The lift is normal to the plane generated by the blade and relative velocity
                lift_unit = np.cross(drag_unit, blade_unit)
                lift_unit_mag = np.linalg.norm(lift_unit)
                if lift_unit_mag < 1e-6:
                    lift_unit_mag = 1e-6
                lift_unit = lift_unit/lift_unit_mag
                
                vector_nodal_drag = np.outer(nodal_drag[:, k], drag_unit)
                vector_nodal_lift = np.outer(nodal_lift[:, k], lift_unit)

                drag_force += vector_nodal_drag
                lift_force += vector_nodal_lift
                    
            # The total turbine force is the sum of lift and drag effects
        turbine_force = drag_force + lift_force

        for k in range(ndim):
            tf_vec[k::ndim] = turbine_force[:, k]
            # tf_lift_vec[k::ndim] += lift_force[:, k]
            # tf_drag_vec[k::ndim] += drag_force[:, k]

        # Save the output
        tf_vec[np.abs(tf_vec) < 1e-12] = 0.0
        # tf_lift_vec[np.abs(tf_lift_vec) < 1e-12] = 0.0
        # tf_drag_vec[np.abs(tf_drag_vec) < 1e-12] = 0.0

        self.problem.tf.vector()[:] = tf_vec

    # ================================================================

    def UpdateActuatorLineForceOld(self, simTime):
        coords = self.problem.fs.V.tabulate_dof_coordinates()
        coords = np.copy(coords[0::self.problem.dom.dim, :])

        # Set up the turbine geometry
        num_blades = 3
        theta_vec = np.linspace(0, 2.0*np.pi, num_blades+1)
        theta_vec = theta_vec[0:num_blades]

        # print(theta_vec)

        # Lift and drag coefficient (could be an array and you interpolate the value based on R)
        cl = 1.0
        cd = 2.0

        rho = 1

        u_inf = 8

        # Blade length (turbine radius)
        L = 13.5

        # Chord length
        c = L/20

        # Number of blade evaluation sections
        num_actuator_nodes = 50
        rdim = np.linspace(0, L, num_actuator_nodes)
        zdim = 0.0 + np.zeros(num_actuator_nodes)
        xblade = np.vstack((np.zeros(num_actuator_nodes), rdim, zdim))

        # Width of individual blade segment
        w = rdim[1] - rdim[0]

        # Width of Gaussian
        eps = 2.5*c

        tf_vec = np.zeros(np.size(coords))




        RPM = 15.0
        period = 60.0/RPM
        theta_offset = simTime/period*2.0*np.pi

        tip_speed = np.pi*2*L*RPM/60

        blade_vel = np.linspace(0.0, tip_speed, num_actuator_nodes)


        constant_vel_mag = True

        if constant_vel_mag:
            # Lift and drag force (calculated outside loop since cl, cd, u_inf, and c assumed constant)
            lift = 0.5*cl*rho*u_inf**2*c*w
            drag = 0.5*cd*rho*u_inf**2*c*w

            # Save time by moving calculation out of loop
            L1 = lift/(eps**3 * np.pi**1.5)
            D1 = drag/(eps**3 * np.pi**1.5)
        else:
            # Lift and drag force (calculated outside loop since cl, cd, u_inf, and c assumed constant)
            u_inf = u_inf**2 + blade_vel**2
            lift = 0.5*cl*rho*c*w*u_inf
            drag = 0.5*cd*rho*c*w*u_inf

            # Save time by moving calculation out of loop
            L1 = lift/(eps**3 * np.pi**1.5)
            D1 = drag/(eps**3 * np.pi**1.5)



        for theta_0 in theta_vec:
            theta = theta_0 + theta_offset
            
            # Create rotation matrix for this turbine blade
            rotA = np.array([[1, 0, 0],
                             [0, np.cos(theta), -np.sin(theta)],
                             [0, np.sin(theta), np.cos(theta)]])

            # Rotate the entire [x; y; z] matrix using this matrix
            xblade_rotated = np.dot(rotA, xblade)
            xblade_rotated[2, :] += 32.1

            use_vectorized_calculation = True

            if use_vectorized_calculation:
                coordsLinear = np.copy(coords.reshape(-1, 1))

                xblade_rotated_full = np.tile(xblade_rotated, (np.shape(coords)[0], 1))

                dx_full = (coordsLinear - xblade_rotated_full)**2

                dist2 = dx_full[0::self.problem.dom.dim] + \
                dx_full[1::self.problem.dom.dim] + \
                dx_full[2::self.problem.dom.dim]

                total_dist2_lift = np.sum(L1*np.exp(-dist2/eps**2), axis = 1)
                total_dist2_drag = np.sum(D1*np.exp(-dist2/eps**2), axis = 1)

                tf_vec[0::self.problem.dom.dim] += -total_dist2_drag
                tf_vec[1::self.problem.dom.dim] += total_dist2_lift*np.sin(theta)
                tf_vec[2::self.problem.dom.dim] += -total_dist2_lift*np.cos(theta)

            else:
                for k, x in enumerate(coords):
                    # Flip row into column
                    xT = x.reshape(-1, 1)

                    # Subtract this 3x1 point from the 3xN array of rotated turbine segments
                    dx = (xT - xblade_rotated)**2
                    mag = np.sum(dx, axis = 0)

                    # Add up the contribution from each blade segment
                    lift_sum = np.sum(L1*np.exp(-mag/eps**2))
                    drag_sum = np.sum(D1*np.exp(-mag/eps**2))

                    # Store the individual vector components using linear index
                    tf_vec[3*k+0] += -drag_sum
                    tf_vec[3*k+1] += lift_sum*np.sin(theta)
                    tf_vec[3*k+2] += -lift_sum*np.cos(theta)
                
        tf_vec[np.abs(tf_vec) < 1e-12] = 0.0

        self.problem.tf.vector()[:] = tf_vec

    # ================================================================

    def UpdateTurbineForce(self, simTime, turbsPerPlatform):
        coords = self.problem.fs.V.tabulate_dof_coordinates()
        coords = np.copy(coords[0::self.problem.dom.dim, :])

        # Pre allocate all numpy arrays and vectors
        tf_array = np.zeros(np.shape(coords))
        tf_vec = np.zeros(np.size(tf_array))
        xs = np.zeros(np.shape(coords))

        # Radius of the two "arms" measured from the hinge
        rad = 189.0

        if turbsPerPlatform == 1:
            rad = 0.0

        # Angle defined between the two "arms"
        phi = np.pi/3.0

        # Calculate the offset from the hinge to each turbine
        xp_offset = rad*np.cos(phi/2.0)
        yp_offset = rad*np.sin(phi/2.0)

        # delta_yaw = 0.0

        for k in range(self.problem.farm.numturbs):
            # Position of the kth turbune
            xpos = float(self.problem.farm.mx[k])
            ypos = float(self.problem.farm.my[k])
            
            if self.problem.dom.dim == 2:
                x0 = np.array([xpos, ypos])
            else:
                zpos = float(self.problem.farm.mz[k])
                x0 = np.array([xpos, ypos, zpos])

            # Yaw, thickness, radius, and mass of the kth turbine
            # If desired, shift each turbine by a constant amount
            # delta_yaw = np.pi/4.0*np.sin(np.pi*(simTime/1000.0 + k/self.problem.farm.numturbs))
            delta_yaw = 0.0

            yaw = float(self.problem.farm.myaw[k] + delta_yaw)
            W = float(self.problem.farm.W[k]/2.0)
            R = float(self.problem.farm.RD[k]/2.0)
            ma = float(self.problem.farm.ma[k])

            # Create a rotation matrix for this yaw angle
            A_rotation = self.RotationMatrix(yaw)

            # Rotate the turbine after shifting (x0, y0, z0) to the center of the turbine
            xs0 = np.dot(coords - x0, A_rotation)

            for doublet in range(turbsPerPlatform):

                offset = np.zeros(self.problem.dom.dim)
                offset[0] = xp_offset
                offset[1] = yp_offset*(-1)**doublet

                # Offset each turbine from the center of rotation
                xs = xs0 - offset

                # Normal to blades: Create the function that represents the Thickness of the turbine
                T_norm = 1.902701539733748
                T = np.exp(-(xs[:, 0]/W)**10.0)/(T_norm*W)

                # Tangential to blades: Create the function that represents the Disk of the turbine
                D_norm = 2.884512175878827
                if self.problem.dom.dim == 2:
                    D1 = (xs[:, 1]/R)**2.0
                else:
                    D1 = (xs[:, 1]/R)**2.0 + (xs[:, 2]/R)**2.0

                D = np.exp(-D1**5.0)/(D_norm*R**2.0)

                # Create the function that represents the force
                if self.problem.dom.dim == 2:
                    r = xs[:, 1]
                else:
                    r = np.sqrt(xs[:, 1]**2.0 + xs[:, 2]**2.0)

                F = 4.0*0.5*(np.pi*R**2.0)*ma/(1.0 - ma)*(r/R*np.sin(np.pi*r/R) + 0.5) * 1.0/.81831

                u_vec = self.problem.u_k1.vector()[:]
                ux = u_vec[0::self.problem.dom.dim]
                uy = u_vec[1::self.problem.dom.dim]
                uD = ux*np.cos(yaw) + uy*np.sin(yaw)

                tf_array[:, 0] = tf_array[:, 0] + F*T*D*np.cos(yaw)*uD**2.0
                tf_array[:, 1] = tf_array[:, 1] + F*T*D*np.sin(yaw)*uD**2.0


        # Riffle shuffle the array elements into a FEniCS-style vector
        for k in range(self.problem.dom.dim):
            tf_vec[k::self.problem.dom.dim] = tf_array[:, k]

        tf_vec[np.abs(tf_vec) < 1e-50] = 0.0

        # Set the vector elements
        self.problem.tf.vector()[:] = tf_vec

    # ================================================================

    def RotationMatrix(self, yaw):
        cosYaw = np.cos(yaw)
        sinYaw = np.sin(yaw)

        if self.problem.dom.dim == 2:
            A_rotation = np.array([[cosYaw, -sinYaw],
                                   [sinYaw,  cosYaw]])
        else:
            A_rotation = np.array([[cosYaw, -sinYaw, 0.0],
                                   [sinYaw,  cosYaw, 0.0],
                                   [   0.0,     0.0, 1.0]])

        return A_rotation

    # ================================================================


    def modifyInletVelocity(self, simTime, bcu):

        # Define tolerance
        tol = 1e-6

        def left_wall(x, on_boundary):
            return on_boundary and x[0] < self.problem.dom.x_range[0] + tol

        HH_vel = self.problem.bd.HH_vel

        # Get the coordinates using the vector funtion space, V
        coords = self.problem.fs.V.tabulate_dof_coordinates()
        coords = np.copy(coords[0::self.problem.dom.dim, :])

        # Create a function representing to the inlet velocity
        vel_inlet_func = Function(self.problem.fs.V)

        inlet_type = 1

        if inlet_type == 1:
            # Create arrays for the steady, vortex, and combined velocities
            vel_steady = np.zeros(np.shape(coords))
            vel_steady[:, 0] = HH_vel
            # print(HH_vel)

            vel_vort = np.zeros(np.shape(coords))
            vel_inlet = np.zeros(np.shape(coords))

            # Specify the vortex radius
            vortRad = 1000
            vortRad2 = vortRad**2

            # Specify the vortex velocity and calculate its position from the starting point
            vortVel = 1.0

            period = 1000.0
            xd = period/2 - vortVel*(simTime%period)

            fac = 0.1
            sep = 650
            Tau = 1000

            for k, x in enumerate(coords):
                if x[0] < self.problem.dom.x_range[0] + tol:

                    # xd should replace x[0] in the following equations
                    if np.abs(xd) < 1e-3:
                        xd = 1e-3

                    cp = ((x[1] + sep/2)**2 + xd**2)/(4*Tau)
                    cn = ((x[1] - sep/2)**2 + xd**2)/(4*Tau)

                    # U-velocity
                    vel_inlet[k, 0] = fac*((1 - np.exp(-cp))/cp*(x[1] + sep/2) -\
                                           (1 - np.exp(-cn))/cn*(x[1] - sep/2)) + 1

                    # V-velocity
                    vel_inlet[k, 1] = fac*(-(1 - np.exp(-cp))/cp*xd +\
                                            (1 - np.exp(-cn))/cn*xd)

                    norm = np.sqrt(vel_inlet[k, 0]*vel_inlet[k, 0] + vel_inlet[k, 1]*vel_inlet[k, 1])

                    if norm > 10.0:
                        vel_inlet[k, 0] = vel_inlet[k, 0]/norm*10.0
                        vel_inlet[k, 1] = vel_inlet[k, 1]/norm*10.0

                    # dx = x - vortPos
                    # dist2 = dx[0]*dx[0] + dx[1]*dx[1]

                    # if dist2 < vortRad2:
                    #     theta = np.arctan2(dx[1], dx[0])
                    #     fac = 1.0 - np.sqrt(dist2/vortRad2)
                    #     vel_vort[k, 0] = -np.sin(theta)
                    #     vel_vort[k, 1] = np.cos(theta)
                    # else:
                    #     fac = 0.0
                    #     vel_vort[k, 0] = 0.0
                    #     vel_vort[k, 1] = 0.0

                    # vel_inlet[k, :] = (1.0-fac)*vel_steady[k, :] + HH_vel*fac*vel_vort[k, :]

        elif inlet_type == 2:
            jet_rad = 400

            vel_inlet = np.zeros(np.shape(coords))

            for k, x in enumerate(coords):
                if x[0] < self.problem.dom.x_range[0] + tol:
                    if np.abs(x[1]) < jet_rad:
                        thetaMax = 15.0/180.0*np.pi

                        theta = thetaMax*np.sin(simTime/1000*2*np.pi)

                        vel_inlet[k, 0] = 2.0*HH_vel*np.cos(theta)
                        vel_inlet[k, 1] = 2.0*HH_vel*np.sin(theta)
                    else:
                        vel_inlet[k, 0] = HH_vel
                        vel_inlet[k, 1] = 0.0


        # Riffle shuffle the array elements into a 1D vector
        vel_inlet_vector = np.zeros(np.size(vel_inlet))

        for k in range(self.problem.dom.dim):
            vel_inlet_vector[k::self.problem.dom.dim] = vel_inlet[:, k]

        # Assign the function the vector of values
        vel_inlet_func.vector()[:] = vel_inlet_vector


        # Update the inlet velocity
        bcu[0] = DirichletBC(self.problem.fs.V, vel_inlet_func, left_wall)

        return bcu


# ================================================================

class MultiAngleSolver(SteadySolver):
    """
    This solver will solve the problem using the steady state solver for every
    angle in angles.

    Args: 
        problem (:meth:`windse.ProblemManager.GenericProblem`): a windse problem object.
        angles (list): A list of wind inflow directions.
    """ 

    def __init__(self,problem):
        super(MultiAngleSolver, self).__init__(problem)
        if self.params["domain"]["type"] in ["imported"]:
            raise ValueError("Cannot use a Multi-Angle Solver with an "+self.params["domain"]["type"]+" domain.")
        self.orignal_solve = super(MultiAngleSolver, self).Solve
        if self.problem.dom.raw_inflow_angle is None:
            self.wind_range = [0, 2.0*np.pi,self.num_wind_angles]
        elif isinstance(self.problem.dom.raw_inflow_angle,list):
            if len(self.problem.dom.raw_inflow_angle)==3:
                self.wind_range = self.problem.dom.raw_inflow_angle
            else:
                self.wind_range = [self.problem.dom.raw_inflow_angle[0],self.problem.dom.raw_inflow_angle[1],self.num_wind_angles]
        else:
            self.wind_range = [self.problem.dom.raw_inflow_angle,self.problem.dom.raw_inflow_angle+2.0*np.pi,self.num_wind_angles]

        self.angles = np.linspace(*self.wind_range,endpoint=self.endpoint)
        self.angles = meteor_to_math(self.angles)
        # self.angles += self.angle_offset

    def Solve(self):
        for i, theta in enumerate(self.angles):
            self.fprint("Performing Solve {:d} of {:d}".format(i+1,len(self.angles)),special="header")
            self.fprint("Wind Angle: "+repr(math_to_meteor(theta)))
            if i > 0 or not near(theta,self.problem.dom.inflow_angle):
                self.problem.dom.inflow_angle = theta
                self.ChangeWindAngle(theta)
            self.iter_val = math_to_meteor(theta)
            self.orignal_solve()
            self.fprint("Finished Solve {:d} of {:d}".format(i+1,len(self.angles)),special="footer")

class TimeSeriesSolver(SteadySolver):
    """
    This solver will solve the problem using the steady state solver for every
    angle in angles.

    Args: 
        problem (:meth:`windse.ProblemManager.GenericProblem`): a windse problem object.
        angles (list): A list of wind inflow directions.
    """ 

    def __init__(self,problem):
        super(TimeSeriesSolver, self).__init__(problem)
        if self.params["domain"]["type"] in ["imported"]:
            raise ValueError("Cannot use a Multi-Angle Solver with an "+self.params["domain"]["type"]+" domain.")
        self.orignal_solve = super(TimeSeriesSolver, self).Solve

        raw_data = np.loadtxt(self.velocity_path,comments="#")
        self.times = raw_data[:,0]
        self.speeds = raw_data[:,1]
        self.angles = raw_data[:,2]
        self.num_solve = len(self.speeds)

    def Solve(self):
        for i in range(self.num_solve):
            time  = self.times[i]
            theta = self.angles[i]
            speed = self.speeds[i]
            self.fprint("Performing Solve {:d} of {:d}".format(i+1,len(self.angles)),special="header")
            self.fprint("Time: "+repr(time))
            self.fprint("Wind Angle: "+repr(theta))
            self.fprint("Wind Speed: "+repr(speed))
            if i > 0 or not near(speed,self.problem.bd.HH_vel):
                self.ChangeWindSpeed(speed)
            if i > 0 or not near(theta,self.problem.dom.inflow_angle):
                self.ChangeWindAngle(theta)
            self.iter_val = time
            self.orignal_solve()

            self.fprint("Finished Solve {:d} of {:d}".format(i+1,len(self.angles)),special="footer")

























##### OLD CODE FOR FUTURE REFERENCE ####


    # # def CalculatePowerFunctional_n(self,inflow_angle = 0.0):
    # def CalculatePowerFunctional(self,inflow_angle = 0.0):

    #     # if self.problem.dom.dim == 3:
    #     if self.problem.dom.dim <= 3:
    #         # print("here")
    #         ### Reconstruct Turbine Force ###
    #         tf = self.problem.tf1*self.u_next[0]**2+self.problem.tf2*self.u_next[1]**2+self.problem.tf3*self.u_next[0]*self.u_next[1]


    #         W = self.problem.farm.thickness[0]*1.0
    #         R = self.problem.farm.RD[0]/2.0 

    #         T_norm = 2.0*gamma(7.0/6.0)
    #         D_norm = pi*gamma(4.0/3.0)
    #         volNormalization_3D = T_norm*D_norm*W*R**2.0
    #         print("3D Volume: "+repr(volNormalization_3D))

    #         T_norm = 2.0*gamma(7.0/6.0)
    #         D_norm = 2.0*gamma(7.0/6.0)
    #         volNormalization_2D = T_norm*D_norm*W*R
    #         print("3D Volume: "+repr(volNormalization_2D))


    #         if self.problem.dom.dim == 2:
    #             dim_scale = volNormalization_3D/volNormalization_2D
    #             print(dim_scale)
    #             dim_scale = pi*R/2
    #             print(dim_scale)
    #             dim_scale = 2*R*R*volNormalization_2D/volNormalization_3D
    #             print(dim_scale)
    #         else: 
    #             dim_scale = 1.0

    #         ### Calculate Power ###
    #         J = dot(tf*dim_scale,self.u_next)*dx

    #         ### Save Individual Powers ###
    #         if self.save_power:
    #             J_list=np.zeros(self.problem.farm.numturbs+1)
    #             J_list[-1]=assemble(J,**self.extra_kwarg)

    #             if self.problem.farm.actuator_disks_list is not None:
    #                 for i in range(self.problem.farm.numturbs):
    #                     yaw = self.problem.farm.myaw[i]+inflow_angle
    #                     tf1 = self.problem.farm.actuator_disks_list[i] * cos(yaw)**2
    #                     tf2 = self.problem.farm.actuator_disks_list[i] * sin(yaw)**2
    #                     tf3 = self.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
    #                     tf = tf1*self.u_next[0]**2+tf2*self.u_next[1]**2+tf3*self.u_next[0]*self.u_next[1]
    #                     J_list[i] = assemble(dot(tf*dim_scale,self.u_next)*dx,**self.extra_kwarg)

    #             for j in J_list:
    #                 print(j)
    #             # exit()
    #     # else:
    #     #     J=0
    #     #     J_list=np.zeros(self.problem.farm.numturbs+1)
    #     #     for i in range(self.problem.farm.numturbs):

    #     #         ### Set some Values ###
    #     #         yaw = self.problem.farm.myaw[i]+inflow_angle
    #     #         W = self.problem.farm.W[i]
    #     #         R = self.problem.farm.RD[i]/2.0 
    #     #         A = pi*R**2.0
    #     #         a = self.problem.farm.ma[i]
    #     #         C_tprime = 4*a/(1-a)
    #     #         C_pprime = 0.45/(1-a)**3

    #     #         ### Define Integral Constants ###
    #     #         T_norm = 2.0*gamma(7.0/6.0)
    #     #         D_norm = 2.0*gamma(7.0/6.0)
    #     #         h1 = float(hyper([],[1/3, 2/3, 5/6, 7/6],-(np.pi**6/46656)))
    #     #         h2 = float(hyper([],[2/3, 7/6, 4/3, 3/2],-(np.pi**6/46656)))
    #     #         h3 = float(hyper([],[4/3, 3/2, 5/3, 11/6],-(np.pi**6/46656)))
    #     #         SD_norm = (1/3*np.pi**(3/2)*h1 - 1/15*np.pi**3*gamma(11/6)*h2 + gamma(7/6)*(1 + 1/360*np.pi**5*h3))/(.81831)

    #     #         ### Reconstruct Turbine Force ###
    #     #         tf1 = self.problem.farm.actuator_disks_list[i] * cos(yaw)**2
    #     #         tf2 = self.problem.farm.actuator_disks_list[i] * sin(yaw)**2
    #     #         tf3 = self.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
    #     #         tf = tf1*self.u_next[0]**2+tf2*self.u_next[1]**2+tf3*self.u_next[0]*self.u_next[1]
    #     #         tf = tf*T_norm*D_norm*W*R

    #     #         ### Calculate Volume of the Domain ###
    #     #         unit = Function(self.problem.fs.Q)
    #     #         unit.vector()[:] = 1.0
    #     #         vol = assemble(unit*dx)

    #     #         ### Calculate Integral of Actuator Disk ###
    #     #         if self.problem.farm.force == "constant":
    #     #             Va = C_tprime*T_norm*D_norm*W*R*R
    #     #         elif self.problem.farm.force == "sine":
    #     #             Va = C_tprime*T_norm*SD_norm*W*R*R

    #     #         ### Calculate Power ###
    #     #         # J_current = 0.5*A*C_pprime*(assemble(F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
    #     #         # J_current = 0.5*A*C_pprime/vol*(dot(tf,self.u_next)/Va)**3*dx
    #     #         J_current = 0.5*A*C_pprime/vol*(dot(tf,self.u_next)/Va)*dx
    #     #         # J_current = 0.5*A*C_pprime/vol*(unit/Va)**3*dx
    #     #         # J_current = unit*0.5*A*C_pprime/vol*(1/Va)**3*dx
    #     #         # J_current = unit*Va*dx
    #     #         J += J_current

    #     #         if self.save_power:
    #     #             J_list[i] = assemble(J_current)

    #     #     if self.save_power:
    #     #         J_list[-1]=np.sum(J_list[:-1])
    #     #         print()
    #     #         for j in J_list:
    #     #             print(j)
    #     #         # exit()

    #         folder_string = self.params.folder+"/data/"
    #         if not os.path.exists(folder_string): os.makedirs(folder_string)

    #         if self.J_saved:
    #             f = open(folder_string+"power_data.txt",'ab')
    #         else:
    #             f = open(folder_string+"power_data.txt",'wb')
    #             self.J_saved = True

    #         np.savetxt(f,[J_list])
    #         f.close()

    #     return J
        
    # # def CalculatePowerFunctional_o(self,delta_yaw = 0.0):
    # # # def CalculatePowerFunctional(self,delta_yaw = 0.0):
    # #     self.fprint("Computing Power Functional")

    # #     x=SpatialCoordinate(self.problem.dom.mesh)
    # #     J=0.
    # #     J_list=np.zeros(self.problem.farm.numturbs+1)
    # #     for i in range(self.problem.farm.numturbs):

    # #         mx = self.problem.farm.mx[i]
    # #         my = self.problem.farm.my[i]
    # #         mz = self.problem.farm.mz[i]
    # #         x0 = [mx,my,mz]
    # #         W = self.problem.farm.thickness[i]*1.0
    # #         R = self.problem.farm.RD[i]/2.0 
    # #         ma = self.problem.farm.ma[i]
    # #         yaw = self.problem.farm.myaw[i]+delta_yaw
    # #         u = self.u_next
    # #         A = pi*R**2.0
    # #         C_tprime = 4*ma/(1-ma)
    # #         C_pprime = 0.45/(1-ma)**3
            
    # #         ### Rotate and Shift the Turbine ###
    # #         xs = self.problem.farm.YawTurbine(x,x0,yaw)
    # #         u_d = u[0]*cos(yaw) + u[1]*sin(yaw)

    # #         ### Create the function that represents the Thickness of the turbine ###
    # #         T = exp(-pow((xs[0]/W),6.0))

    # #         if self.problem.dom.dim == 3:
    # #             # WTGbase = Expression(("cos(yaw)","sin(yaw)","0.0"),yaw=yaw,degree=1)
    # #             WTGbase = as_vector((cos(yaw),sin(yaw),0.0))

    # #             ### Create the function that represents the Disk of the turbine
    # #             r = sqrt(xs[1]**2.0+xs[2]**2.0)/R
    # #             D = exp(-pow(r,6.0))
    # #             # D = exp(-pow((pow((xs[1]/R),2)+pow((xs[2]/R),2)),6.0))

    # #             volNormalization = assemble(T*D*dx)
    # #             print(volNormalization)
    # #             T_norm = 2.0*gamma(7.0/6.0)
    # #             D_norm = pi*gamma(4.0/3.0)
    # #             volNormalization = T_norm*D_norm*W*R**2.0
    # #             print(volNormalization)

    # #             ### Create the function that represents the force ###
    # #             if self.problem.farm.force == "constant":
    # #                 F = 0.5*A*C_tprime
    # #             elif self.problem.farm.force == "sine":
    # #                 # r = sqrt(xs[1]**2.0+xs[2]**2)
    # #                 F = 0.5*A*C_tprime*(r*sin(pi*r)+0.5)/(.81831)

    # #             J_current = dot(F*T*D*WTGbase*u_d**2,u)*dx
    # #             J += J_current
    # #             if self.save_power:
    # #                 J_list[i] = assemble(J_current)

    # #         else:
    # #             # WTGbase = Expression(("cos(yaw)","sin(yaw)"),yaw=yaw,degree=1)
    # #             WTGbase = as_vector((cos(yaw),sin(yaw)))

    # #             ### Create the function that represents the Disk of the turbine
    # #             r = sqrt(xs[1]**2.0+xs[2]**2.0)/R
    # #             D = exp(-pow(r,6.0))


    # #             # T_norm = 2.0*gamma(7.0/6.0)
    # #             # D_norm = 2.0*gamma(7.0/6.0)
    # #             # h1 = float(hyper([],[1/3, 2/3, 5/6, 7/6],-(np.pi**6/46656)))
    # #             # h2 = float(hyper([],[2/3, 7/6, 4/3, 3/2],-(np.pi**6/46656)))
    # #             # h3 = float(hyper([],[4/3, 3/2, 5/3, 11/6],-(np.pi**6/46656)))
    # #             # SD_norm = (1/3*np.pi**(3/2)*h1 - 1/15*np.pi**3*gamma(11/6)*h2 + gamma(7/6)*(1 + 1/360*np.pi**5*h3))/(.81831)
    # #             # volNormalization = T_norm*SD_norm*W*R**(self.problem.dom.dim)
    # #             volNormalization = assemble(T*D*dx)

    # #             T_norm = 2.0*gamma(7.0/6.0)
    # #             D_norm = pi*gamma(4.0/3.0)
    # #             volNormalization_3D = T_norm*D_norm*W*R**2.0
    # #             # print(volNormalization_3D)


    # #             T_norm = 2.0*gamma(7.0/6.0)
    # #             D_norm = 2.0*gamma(7.0/6.0)
    # #             volNormalization_2D = T_norm*D_norm*W*R
    # #             # print(volNormalization_2D)


    # #             print(volNormalization,volNormalization_2D,volNormalization_3D,2*volNormalization_2D/volNormalization_3D)

    # #             ### Create the function that represents the force ###
    # #             if self.problem.farm.force == "constant":
    # #                 F = 0.5*self.problem.farm.RD[i]*C_tprime    
    # #             elif self.problem.farm.force == "sine":
    # #                 F = 0.5*self.problem.farm.RD[i]*C_tprime*(r*sin(pi*r)+0.5)/(.81831)

    # #             # V  = assemble(F*T*D*dx)
    # #             # Va = float(C_tprime)*T_norm*SD_norm*W*R*R
    # #             # print(V,Va,V/Va,V/(W*R*R),Va/(W*R*R))

    # #             J_current = assemble(dot(F*T*D*WTGbase*u_d**2,u)*dx)
    # #             # J_current = (assemble(dot(F*T*D*WTGbase*u_d**2,u)*dx)
    # #             # J_current = 0.5*A*C_pprime*(assemble(dot(F*T*D*WTGbase*u_d**2,u)*dx)/assemble(F*T*D*dx))
    # #             # print(float(J_current))
    # #             # J_current_old = 0.5*A*C_pprime*(assemble(F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
    # #             # J_current = 0.5*A*C_pprime*(assemble(T*D*u*dx)/assemble(T*D*dx))**3

    # #             # J_3D = [1692363.167076575,801778.7751333286,545135.3528982735]
    # #             # J_3D = [767698.3159772983,420644.4831798226,291267.477222329]
    # #             # print(float(J_current_old),float(J_current_old/J_current),float(J_3D[i]/J_current))
    # #             # J_current = (assemble(0.5*A*C_pprime*F*T*D*u_d*dx)/assemble(F*T*D*dx))




    # #             ### LIST O PROBLEMS ###
    # #             # 1. Unfortunately, this is not a form so cannot be used with dolfin_adjoint (this just returns a float)
    # #             # 2. Where did WTGbase go? I know it need to be a scaler but does yaw angle not matter in 2D?
    # #             # 3. Should the C_pprime term be inside the numerator? 
    # #             # 4. Are you positive you want to divide by the total force (F*T*D) instead of just the space kernal (T*D)



 










    # #             # J_current = 0.5*A*C_pprime*(1/assemble(F*T*D*dx))**3
    # #             J += J_current
    # #             if self.save_power:
    # #                 # J_list[i] = (assemble(F*T*D*u_d*dx)/assemble(F*T*D*dx))**3
    # #                 J_list[i] = J_current

    # #     if self.save_power:
    # #         J_list[-1]=np.sum(J_list[:-1])
    # #         print()
    # #         for j in J_list:
    # #             print(j)
    # #         exit()

    # #         folder_string = self.params.folder+"/data/"
    # #         if not os.path.exists(folder_string): os.makedirs(folder_string)

    # #         if self.J_saved:
    # #             f = open(folder_string+"power_data.txt",'ab')
    # #         else:
    # #             f = open(folder_string+"power_data.txt",'wb')
    # #             self.J_saved = True

    # #         np.savetxt(f,[J_list])
    # #         f.close()

    # #     return J

