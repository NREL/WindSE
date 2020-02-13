"""
The SolverManager contains all the different ways to solve problems generated
in windse
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    from sys import platform
    import time
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
    # from memory_profiler import memory_usage

    ### Import the cumulative parameters ###
    from windse import windse_parameters, CalculateActuatorLineTurbineForces

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

        ### Import objective functions ###
        import windse.ObjectiveFunctions as obj_funcs

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt

    ### Improve Solver parameters ###
    parameters["std_out_all_processes"] = False;
    parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -fno-math-errno -march=native'        
    parameters["form_compiler"]["optimize"]     = True
    parameters["form_compiler"]["cpp_optimize"] = True
    # parameters['form_compiler']['representation'] = 'tsfc'
    parameters['form_compiler']['representation'] = 'uflacs'
    if windse_parameters["wind_farm"].get("turbine_space","Quadrature") == "Quadrature":
        parameters['form_compiler']['quadrature_degree'] = windse_parameters["wind_farm"].get("turbine_degree",6)
    else:
        parameters['form_compiler']['quadrature_degree'] = 6

class GenericSolver(object):
    """
    A GenericSolver contains on the basic functions required by all solver objects.
    """
    def __init__(self,problem):
        self.params = windse_parameters
        self.problem  = problem
        # self.u_next,self.p_next = self.problem.up_next.split(True)
        self.u_next,self.p_next = split(self.problem.up_next)
        self.nu_T = self.problem.nu_T
        self.first_save = True
        self.fprint = self.params.fprint
        self.extra_kwarg = {}
        if self.params["general"].get("dolfin_adjoint", False):
            self.extra_kwarg["annotate"] = False
        self.save_power = self.params["solver"].get("save_power",False)

        #Check if we are optimizing
        if self.params.get("optimization",{}):
            self.optimizing = True
            self.J = 0
            self.objective_type = self.params["optimization"]["objective_type"]
            self.objective_func = obj_funcs.objectives_dict[self.objective_type]
        else:
            self.optimizing = False

        #Check if we need to save the power output
        if self.save_power:
            self.J = 0
            self.J_saved = False

    # def Plot(self):
    #     """
    #     This function plots the solution functions using matplotlib and saves the 
    #     output to output/.../plots/u.pdf and output/.../plots/p.pdf
    #     """

    #     ### Create the path names ###
    #     folder_string = self.params.folder+"/plots/"
    #     u_string = self.params.folder+"/plots/u.pdf"
    #     p_string = self.params.folder+"/plots/p.pdf"

    #     ### Check if folder exists ###
    #     if not os.path.exists(folder_string): os.makedirs(folder_string)

    #     ### Plot the x component of velocity ###
    #     plot(self.u_next[0],title="Velocity in the x Direction")
    #     plt.savefig(u_string)
    #     plt.figure()

    #     ### Plot the pressure ###
    #     plot(self.p_next,title="Pressure")
    #     plt.savefig(p_string)
    #     plt.show()

    def Save(self,val=0):
        """
        This function saves the mesh and boundary markers to output/.../solutions/
        """
        u,p = self.problem.up_next.split(True,**self.extra_kwarg)
        u.vector()[:]=u.vector()[:]/self.problem.dom.xscale
        self.problem.dom.mesh.coordinates()[:]=self.problem.dom.mesh.coordinates()[:]/self.problem.dom.xscale

        if self.first_save:
            self.u_file = self.params.Save(u,"velocity",subfolder="solutions/",val=val)
            self.p_file = self.params.Save(p,"pressure",subfolder="solutions/",val=val)
            # self.nuT_file = self.params.Save(self.nu_T,"eddy_viscosity",subfolder="solutions/",val=val)
            self.first_save = False
        else:
            self.params.Save(u,"velocity",subfolder="solutions/",val=val,file=self.u_file)
            self.params.Save(p,"pressure",subfolder="solutions/",val=val,file=self.p_file)
            # self.params.Save(self.nu_T,"eddy_viscosity",subfolder="solutions/",val=val,file=self.nuT_file)
        u.vector()[:]=u.vector()[:]*self.problem.dom.xscale
        self.problem.dom.mesh.coordinates()[:]=self.problem.dom.mesh.coordinates()[:]*self.problem.dom.xscale

    def ChangeWindSpeed(self,speed):
        """
        This function recomputes all necessary components for a new wind direction

        Args: 
            theta (float): The new wind angle in radians
        """
        self.problem.ChangeWindSpeed(speed)

    def ChangeWindAngle(self,theta):
        """
        This function recomputes all necessary components for a new wind direction

        Args: 
            theta (float): The new wind angle in radians
        """
        self.problem.ChangeWindAngle(theta)


    def CalculatePowerFunctional(self,inflow_angle = 0.0):
        J = -assemble(dot(self.problem.tf,self.u_next)*dx)
        return J

    def SavePower(self,inflow_angle=0.0):

        J_list=np.zeros(self.problem.farm.numturbs+1)

        if self.problem.farm.actuator_disks_list is not None:
            for i in range(self.problem.farm.numturbs):
                yaw = self.problem.farm.myaw[i]+inflow_angle
                tf1 = self.problem.farm.actuator_disks_list[i] * cos(yaw)**2
                tf2 = self.problem.farm.actuator_disks_list[i] * sin(yaw)**2
                tf3 = self.problem.farm.actuator_disks_list[i] * 2.0 * cos(yaw) * sin(yaw)
                tf = tf1*self.u_next[0]**2+tf2*self.u_next[1]**2+tf3*self.u_next[0]*self.u_next[1]
                J_list[i] = assemble(dot(tf,self.u_next)*dx,**self.extra_kwarg)
        J_list[-1]=sum(J_list)

        folder_string = self.params.folder+"/data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if self.J_saved:
            f = open(folder_string+"power_data.txt",'ab')
        else:
            f = open(folder_string+"power_data.txt",'wb')
            self.J_saved = True

        np.savetxt(f,[J_list])
        f.close()

class SteadySolver(GenericSolver):
    """
    This solver is for solving the steady state problem

    Args: 
        problem (:meth:`windse.ProblemManager.GenericProblem`): a windse problem object.
    """
    def __init__(self,problem):
        super(SteadySolver, self).__init__(problem)

    def Solve(self,iter_val=0):
        """
        This solves the problem setup by the problem object.
        """

        ### Save Files before solve ###
        self.fprint("Saving Input Data",special="header")
        if "mesh" in self.params.output:
            self.problem.dom.Save(val=iter_val)
        if "initial_guess" in self.params.output:
            self.problem.bd.SaveInitialGuess(val=iter_val)
        if "height" in self.params.output and self.problem.dom.dim == 3:
            self.problem.bd.SaveHeight(val=iter_val)
        if "turbine_force" in self.params.output:
            self.problem.farm.SaveActuatorDisks(val=iter_val)
        self.fprint("Finished",special="footer")

        ####################################################################
        ### This is the better way to define a nonlinear problem but it 
        ### doesn't play nice with dolfin_adjoint
        ### Define Jacobian ###
        # dU = TrialFunction(self.problem.fs.W)
        # J  = derivative(self.problem.F,  self.problem.up_next, dU)

        # ### Setup nonlinear solver ###
        # nonlinear_problem = NonlinearVariationalProblem(self.problem.F, self.problem.up_next, self.problem.bd.bcs, J)
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

        use_fast_solver = True

        if use_fast_solver:
            # ### Add some helper functions to solver options ###


            krylov_options = {"absolute_tolerance": 9e-3,
                              "relative_tolerance": 9e-1,
                              "maximum_iterations": 5000}

            newton_options = {"relaxation_parameter": 0.9,
                              "maximum_iterations": 20,
                              "linear_solver": "mumps",
                              # "preconditioner": "none",
                              "absolute_tolerance": 1e-6,
                              "relative_tolerance": 1e-5,
                              "krylov_solver": krylov_options}

            solver_parameters = {"nonlinear_solver": "newton",
                                 "newton_solver": newton_options}
                                 # {'ksp_gmres_restart': 100}

        else:
            solver_parameters = {"nonlinear_solver": "snes",
                                 "snes_solver": {
                                 "absolute_tolerance": 1e-6,
                                 "relative_tolerance": 1e-5,
                                 "linear_solver": "mumps", 
                                 "maximum_iterations": 40,
                                 "error_on_nonconvergence": True,
                                 "line_search": "bt", # Should this should be changed to "basic"?
                                 # "line_search": "basic",
                                 }}

        ### Start the Solve Process ###
        self.fprint("Solving",special="header")
        start = time.time()
        
        # ### Solve the Baseline Problem ###
        # solve(self.problem.F_sans_tf == 0, self.problem.up_next, self.problem.bd.bcs, solver_parameters=solver_parameters, **self.extra_kwarg)

        # ### Store the Baseline and Assign for the real solve ###
        # self.up_baseline = self.problem.up_next.copy(deepcopy=True)
        # self.problem.up_next.assign(self.up_baseline)

        ### Solve the real problem ###
        # mem0=memory_usage()[0]
        # mem_out, _ = memory_usage((solve,(self.problem.F == 0, self.problem.up_next, self.problem.bd.bcs),{"solver_parameters": solver_parameters}),max_usage=True,retval=True,max_iterations=1)
        solve(self.problem.F == 0, self.problem.up_next, self.problem.bd.bcs, solver_parameters=solver_parameters)
        stop = time.time()

        self.fprint("Solve Complete: {:1.2f} s".format(stop-start),special="footer")
        # self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))
        # self.u_next,self.p_next = self.problem.up_next.split(True)
        self.u_next,self.p_next = split(self.problem.up_next)
        # self.nu_T = project(self.problem.nu_T,self.problem.fs.Q,solver_type='mumps',**self.extra_kwarg)
        self.nu_T = None

        ### Save solutions ###
        if "solution" in self.params.output:
            self.fprint("Saving Solution",special="header")
            self.Save(val=iter_val)
            self.fprint("Finished",special="footer")

        ### calculate the power for each turbine ###
        ###################################
        ### Fix how angle is transfered ###
        ###################################
        if self.save_power:
            self.SavePower(((iter_val-self.problem.dom.init_wind)))


        if self.optimizing:
            self.J += self.objective_func(self,(iter_val-self.problem.dom.init_wind)) 
            print(self.outflow_markers)
            # self.J += -dot(self.problem.farm.rotor_disks,self.u_next)*dx

        # self.fprint("Speed Percent of Inflow Speed")
        # ps = []
        # for i in range(6):
        #     HH = self.problem.farm.HH[0]
        #     RD = self.problem.farm.RD[0]
        #     x_val = (i+1)*RD
        #     vel = self.problem.up_next([x_val,0,HH])
        #     vel = vel[0:3]
        #     nom = np.linalg.norm(vel)
        #     perc = nom/self.problem.bd.HH_vel
        #     ps.append(perc)
        #     self.fprint("Speed Percent at ("+repr(int(x_val))+", 0, "+repr(HH)+"): "+repr(perc))
        # print(ps)

# 
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

    def Solve(self,iter_val=0):
        # Start the unsteady solver ONLY if an unsteady problem has been created
        if self.problem.params["problem"]["type"] == 'unsteady':
            self.fprint("Solving with UnsteadySolver", special="header")
        else:
            raise ValueError("UnsteadySolver can only be run with ProblemType = unsteady, not %s" \
                % (self.problem.params["problem"]["type"]))

        # ================================================================

        # Define the final simulation time
        # FIXME: This should also be set in params.yaml input file
        # tFinal = 6000.0
        tFinal = self.params["solver"].get("final_time", 1)
        
        # Start a counter for the total simulation time
        simTime = 0.0
        recordTime = 4.0
        # recordTime = 1.5*(self.problem.dom.x_range[1]/self.problem.bd.HH_vel)
        # if tFinal < recordTime + 60.0/self.problem.rpm:
        #     self.fprint("Warning: Final time is too small... overriding")
        #     tFinal = recordTime + 60.0/self.problem.rpm
        #     self.fprint("         New Final Time: {:1.2f} s".format(tFinal))
        self.record_delta = tFinal-recordTime

        self.fprint("dt: %.4f" % (self.problem.dt))
        self.fprint("tFinal: %.1f" % (tFinal))

        # ================================================================

        # Specify how frequently to save output files
        saveInterval = self.params["solver"].get("save_interval",1)

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
        self.SaveTimeSeries(simTime)

        self.fprint("Saving Input Data",special="header")
        if "mesh" in self.params.output:
            self.problem.dom.Save(val=iter_val)
        if "initial_guess" in self.params.output:
            self.problem.bd.SaveInitialGuess(val=iter_val)
        if "height" in self.params.output and self.problem.dom.dim == 3:
            self.problem.bd.SaveHeight()
        # if "turbine_force" in self.params.output:
        #     self.problem.farm.SaveTurbineForce(val=iter_val)

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

        # ================================================================

        print('Using Actuator/Turbine Parameters...')
        print('Hub Height: ', self.problem.farm.HH[0])
        print('Yaw: ', self.problem.farm.yaw[0])
        print('Radius: ', self.problem.farm.radius[0])

        self.fprint("Solving",special="header")
        self.fprint("Sim Time | Next dt | U_max")
        self.fprint("--------------------------")

        start = time.time()

        dfd_c_lift = CalculateActuatorLineTurbineForces(self.problem, simTime, dfd='c_lift')
        print('dfd_c_lift: ', np.shape(dfd_c_lift))
        dfd_c_drag = CalculateActuatorLineTurbineForces(self.problem, simTime, dfd='c_drag')
        print('dfd_c_drag: ', np.shape(dfd_c_drag))
        coords = self.problem.fs.V.tabulate_dof_coordinates()
        coords = np.copy(coords[0::self.problem.dom.dim, :])
        print(np.shape(coords))

        exit()

        while simTime < tFinal:
            # Get boundary conditions specific to this timestep
            # bcu, bcp = self.GetBoundaryConditions(simTime/tFinal)
            # bcu = self.modifyInletVelocity(simTime, bcu)

            
            # self.UpdateActuatorLineForce(simTime) # Single turbine, line actuator

            # self.problem.tf = self.problem.farm.TurbineForce_numpy(None,None,None)
            # self.UpdateTurbineForce(simTime, 1) # Single turbine, disk actuator
            # self.UpdateTurbineForce(simTime, 2) # Dubs
            # t1 = time.time()
            # t2 = time.time()
            # print(t2-t1)

            self.problem.bd.UpdateVelocity(simTime)

            # Record the "old" max velocity (before this update)
            u_max_k1 = self.problem.u_k.vector().max()

            # Step 1: Tentative velocity step
            b1 = assemble(self.problem.L1, tensor=b1)
            [bc.apply(b1) for bc in self.problem.bd.bcu]
            solve(A1, self.problem.u_k.vector(), b1, 'gmres', 'default')

            # Step 2: Pressure correction step
            b2 = assemble(self.problem.L2, tensor=b2)
            [bc.apply(b2) for bc in self.problem.bd.bcp]
            solve(A2, self.problem.p_k.vector(), b2, 'gmres', 'hypre_amg')

            # Step 3: Velocity correction step
            b3 = assemble(self.problem.L3, tensor=b3)
            solve(A3, self.problem.u_k.vector(), b3, 'gmres', 'default')

            # Old <- New update step
            self.problem.u_k2.assign(self.problem.u_k1)
            self.problem.u_k1.assign(self.problem.u_k)
            self.problem.p_k1.assign(self.problem.p_k)

            # Record the updated max velocity
            u_max = self.problem.u_k.vector().max()

            # Update the simulation time
            simTime += self.problem.dt

            # Update the turbine force
            new_tf = CalculateActuatorLineTurbineForces(self.problem, simTime)
            self.problem.tf.assign(new_tf)

            # Calculate the objective function
            if self.optimizing and simTime >= recordTime:
                self.J += self.problem.dt/(tFinal-recordTime)*self.objective_func(self,(iter_val-self.problem.dom.init_wind)) 


            if save_next_timestep:
                # Read in new inlet values
                # bcu = self.updateInletVelocityFromFile(saveCount, bcu)
                
                # Clean up simTime to avoid accumulating round-off error
                saveCount += 1
                simTime = saveInterval*saveCount

                # Save output files
                # self.SaveTimeSeries(fp, simTime)
                self.SaveTimeSeries(simTime)

            # Adjust the timestep size, dt, for a balance of simulation speed and stability
            save_next_timestep = self.AdjustTimestepSize(save_next_timestep, saveInterval, simTime, u_max, u_max_k1)

            # After changing timestep size, A1 must be reassembled
            # FIXME: This may be unnecessary (or could be sped up by changing only the minimum amount necessary)
            A1 = assemble(self.problem.a1, tensor=A1)
            [bc.apply(A1) for bc in self.problem.bd.bcu]

            # Print some solver statistics
            self.fprint("%8.2f | %7.2f | %5.2f" % (simTime, self.problem.dt, u_max))

        stop = time.time()

        self.fprint("Finished",special="footer")
        self.fprint("Solve Complete: {:1.2f} s".format(stop-start),special="footer")

    # ================================================================

    def SaveTimeSeries(self, simTime):

        if self.first_save:
            self.velocity_file = self.params.Save(self.problem.u_k,"velocity",subfolder="timeSeries/",val=simTime)
            self.pressure_file   = self.params.Save(self.problem.p_k,"pressure",subfolder="timeSeries/",val=simTime)
            self.turb_force_file   = self.params.Save(self.problem.tf,"turbine_force",subfolder="timeSeries/",val=simTime)
            self.first_save = False
        else:
            self.params.Save(self.problem.u_k,"velocity",subfolder="timeSeries/",val=simTime,file=self.velocity_file)
            self.params.Save(self.problem.p_k,"pressure",subfolder="timeSeries/",val=simTime,file=self.pressure_file)
            self.params.Save(self.problem.tf,"turbine_force",subfolder="timeSeries/",val=simTime,file=self.turb_force_file)

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

        # Set the CFL target (0.2 is a good value for stability and speed, YMMV)
        cfl_target = 0.2

        # Enforce a minimum timestep size
        dt_min = 0.01

        # Calculate the change in velocity using a first-order, backward difference
        dudt = u_max - u_max_k1

        # Calculate the projected velocity
        u_max_projected = u_max + dudt

        # Calculate the ideal timestep size (ignore file output considerations for now)
        dt_new = cfl_target * self.problem.dom.mesh.hmin() / u_max_projected

        # Move to larger dt slowly (smaller dt happens instantly)
        if dt_new > self.problem.dt:
            # Amount of new dt to use: 0 = none, 1 = all
            SOR = 0.5
            dt_new = SOR*dt_new + (1.0-SOR)*self.problem.dt

        # Calculate the time remaining until the next file output
        time_remaining = saveInterval - (simTime % saveInterval)

        # If the new timestep would jump past a save point, modify the new timestep size
        if not save_next_timestep and dt_new + dt_min >= time_remaining:
            dt_new = time_remaining
            save_next_timestep = True
        else:
            save_next_timestep = False

        # dt_new = 0.04

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
        num_blade_segments = 10

        #================================================================
        # Set Derived Constants
        #================================================================

        # Calculate the blade velocity
        period = 60.0/RPM
        tip_speed = np.pi*2.0*L*RPM/60.0
        blade_vel = np.vstack((np.zeros(num_blade_segments),
                               np.zeros(num_blade_segments),
                               np.linspace(0.0, tip_speed, num_blade_segments)))

        # Set the initial angle of each blade
        theta_vec = np.linspace(0.0, 2.0*np.pi, num_blades+1)
        theta_vec = theta_vec[0:num_blades]

        # Calculate discrete node positions
        rdim = np.linspace(0.0, L, num_blade_segments)

        # Calculate width of individual blade segment
        w = rdim[1] - rdim[0]

        # Calculate an array describing the x, y, z position of each point
        xblade = np.vstack((np.zeros(num_blade_segments),
                            rdim,
                            np.zeros(num_blade_segments)))

        #================================================================
        # Begin Calculating Turbine Forces
        #================================================================

        # Lift and drag coefficient (could be an array and you interpolate the value based on R)
        # cl_dolf = Constant((np.linspace(1.5, 0.5, num_blade_segments)))
        # cd_dolf = Constant((np.ones(num_blade_segments)))
        # cl = cl_dolf.values()
        # cd = cd_dolf.values()

        cl = np.linspace(0.0, 2.0, num_blade_segments) # Uncomment for controllability study
        cd = np.linspace(2.0, 0.0, num_blade_segments)
        # cl = np.linspace(2.0, 0.0, num_blade_segments) # Uncomment for controllability study
        # cd = np.linspace(0.0, 2.0, num_blade_segments)
        # cl = np.ones(num_blade_segments)
        # cd = np.ones(num_blade_segments)


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

            # Tile the blade coordinates for every mesh point, [numGridPts*ndim x num_blade_segments]
            xblade_rotated_full = np.tile(xblade_rotated, (np.shape(coords)[0], 1))

            # Subtract and square to get the dx^2 values in the x, y, and z directions
            dx_full = (coordsLinear - xblade_rotated_full)**2

            # Add together to get |x^2 + y^2 + z^2|^2
            dist2 = dx_full[0::ndim] + dx_full[1::ndim] + dx_full[2::ndim]

            # Set if using local velocity around inidividual nodes
            using_local_velocity = False
        
            if using_local_velocity:
                # Generate the fluid velocity from the actual node locations in the flow
                u_fluid = np.zeros((3, num_blade_segments))
                
                for k in range(num_blade_segments):
                    u_fluid[:, k] = self.problem.u_k1(xblade_rotated[0, k],
                                                      xblade_rotated[1, k],
                                                      xblade_rotated[2, k])
                                    
            else:
                # Generate the fluid velocity analytically using the hub height velocity
                # u_inf_vec = u_inf*np.ones(num_blade_segments)
                
                # u_fluid = np.vstack((u_inf_vec,
                #                      np.zeros(num_blade_segments),
                #                      np.zeros(num_blade_segments)))
                u_fluid = np.zeros((3, num_blade_segments))
                
                for k in range(num_blade_segments):
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
            
            for k in range(num_blade_segments):
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
        num_blade_segments = 50
        rdim = np.linspace(0, L, num_blade_segments)
        zdim = 0.0 + np.zeros(num_blade_segments)
        xblade = np.vstack((np.zeros(num_blade_segments), rdim, zdim))

        # Width of individual blade segment
        w = rdim[1] - rdim[0]

        # Width of Gaussian
        eps = 2.5*c

        tf_vec = np.zeros(np.size(coords))




        RPM = 15.0
        period = 60.0/RPM
        theta_offset = simTime/period*2.0*np.pi

        tip_speed = np.pi*2*L*RPM/60

        blade_vel = np.linspace(0.0, tip_speed, num_blade_segments)


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
        self.wind_range = self.params["solver"].get("wind_range", None)
        if  self.wind_range is None:
            self.wind_range = [0, 2.0*np.pi]
            self.endpoint = self.params["solver"].get("endpoint", False)
        else:
            self.endpoint = self.params["solver"].get("endpoint", True)

        self.num_wind = self.params["solver"]["num_wind_angles"]
        self.angle_offset = self.params["solver"].get("angle_offset", 0.0)

        self.angles = np.linspace(self.wind_range[0],self.wind_range[1],self.num_wind,endpoint=self.endpoint)
        self.angles += self.angle_offset

    def Solve(self):
        for i, theta in enumerate(self.angles):
            self.fprint("Performing Solve {:d} of {:d}".format(i+1,len(self.angles)),special="header")
            self.fprint("Wind Angle: "+repr(theta))
            if i > 0 or not near(theta,self.problem.dom.init_wind):
                self.ChangeWindAngle(theta)
            self.orignal_solve(iter_val=theta)
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
        self.velocity_path = self.params["solver"]["velocity_path"]


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
            if i > 0 or not near(theta,self.problem.dom.init_wind):
                self.ChangeWindAngle(theta)
            self.orignal_solve(iter_val=time)

            self.fprint("Finished Solve {:d} of {:d}".format(i+1,len(self.angles)),special="footer")