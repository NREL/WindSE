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

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *

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
    parameters['form_compiler']['representation'] = 'uflacs'
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
        else:
            self.optimizing = False

        #Check if we need to save the power output
        if self.save_power:
            self.J = 0
            self.J_saved = False

    def Plot(self):
        """
        This function plots the solution functions using matplotlib and saves the 
        output to output/.../plots/u.pdf and output/.../plots/p.pdf
        """

        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        u_string = self.params.folder+"/plots/u.pdf"
        p_string = self.params.folder+"/plots/p.pdf"

        ### Check if folder exists ###
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        ### Plot the x component of velocity ###
        plot(self.u_next[0],title="Velocity in the x Direction")
        plt.savefig(u_string)
        plt.figure()

        ### Plot the pressure ###
        plot(self.p_next,title="Pressure")
        plt.savefig(p_string)
        plt.show()

    def Save(self,val=0):
        """
        This function saves the mesh and boundary markers to output/.../solutions/
        """
        u,p = self.problem.up_next.split(True,**self.extra_kwarg)
        if self.first_save:
            self.u_file = self.params.Save(u,"velocity",subfolder="solutions/",val=val)
            self.p_file = self.params.Save(p,"pressure",subfolder="solutions/",val=val)
            # self.nuT_file = self.params.Save(self.nu_T,"eddy_viscosity",subfolder="solutions/",val=val)
            self.first_save = False
        else:
            self.params.Save(u,"velocity",subfolder="solutions/",val=val,file=self.u_file)
            self.params.Save(p,"pressure",subfolder="solutions/",val=val,file=self.p_file)
            # self.params.Save(self.nu_T,"eddy_viscosity",subfolder="solutions/",val=val,file=self.nuT_file)

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
        
    def CalculatePowerFunctional(self,delta_yaw = 0.0):
        self.fprint("Computing Power Functional")

        x=SpatialCoordinate(self.problem.dom.mesh)
        J=0.
        J_list=np.zeros(self.problem.farm.numturbs+1)
        for i in range(self.problem.farm.numturbs):

            mx = self.problem.farm.mx[i]
            my = self.problem.farm.my[i]
            mz = self.problem.farm.mz[i]
            x0 = [mx,my,mz]
            W = self.problem.farm.W[i]*1.0
            R = self.problem.farm.RD[i]/2.0 
            ma = self.problem.farm.ma[i]
            yaw = self.problem.farm.myaw[i]+delta_yaw
            u = self.u_next
            A = pi*R**2.0
            C_tprime = 4*ma/(1-ma)


            ### Rotate and Shift the Turbine ###
            xs = self.problem.farm.YawTurbine(x,x0,yaw)

            if self.problem.dom.dim == 3:
                # WTGbase = Expression(("cos(yaw)","sin(yaw)","0.0"),yaw=yaw,degree=1)
                WTGbase = as_vector((cos(yaw),sin(yaw),0.0))

                ### Create the function that represents the Disk of the turbine
                D_norm = 2.914516237206873
                D = exp(-pow((pow((xs[1]/R),2)+pow((xs[2]/R),2)),6.0))/(D_norm*R**2.0)
            else:
                # WTGbase = Expression(("cos(yaw)","sin(yaw)"),yaw=yaw,degree=1)
                WTGbase = as_vector((cos(yaw),sin(yaw)))
                ### Create the function that represents the Disk of the turbine
                D_norm = 1.916571364345665
                D = exp(-pow((pow((xs[1]/R),2)),6.0))/(D_norm*R**2.0)

            ### Create the function that represents the Thickness of the turbine ###
            T_norm = 1.855438667500383
            T = exp(-pow((xs[0]/W),6.0))/(T_norm*W)
            
            ### Create the function that represents the force ###
            if self.problem.farm.force == "constant":
                F = 0.5*A*C_tprime
            elif self.problem.farm.force == "sine":
                r = sqrt(xs[1]**2.0+xs[2]**2)
                F = 0.5*A*C_tprime*(r/R*sin(pi*r/R)+0.5)/(.81831)
            
            u_d = u[0]*cos(yaw) + u[1]*sin(yaw)

            # J += dot(A*T*D*WTGbase*u_d**2.0,u)*dx
            # J += dot(F*T*D*WTGbase*u_d**2.0,u)*dx
            # J += 0.5*A*C_tprime*((F*T*D*u_d*dx)/(F*T*D*dx))**3



            if self.save_power:
                # J_list[i] = assemble(dot(F*T*D*WTGbase*u_d**2.0,u)*dx)
                J_list[i] = 0.5*A*C_tprime*(assemble(F*T*D*u_d*dx)/assemble(F*T*D*dx))**3

        
        if self.save_power:
            J_list[-1]=np.sum(J_list[:-1])

            folder_string = self.params.folder+"/data/"
            if not os.path.exists(folder_string): os.makedirs(folder_string)

            if self.J_saved:
                f = open(folder_string+"power_data.txt",'ab')
            else:
                f = open(folder_string+"power_data.txt",'wb')
                self.J_saved = True

            np.savetxt(f,[J_list])
            f.close()

        return J_list[-1]

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
            self.problem.farm.SaveRotorDisks(val=iter_val)
        self.fprint("Finished",special="footer")

        ####################################################################
        ### This is the better way to define a nonlinear problem but it 
        ### doesn't play nice with dolfin_adjoint
        # ### Define Jacobian ###
        # dU = TrialFunction(self.problem.fs.W)
        # J  = derivative(self.problem.F,  self.problem.up_next, dU)

        # ### Setup nonlinear solver ###
        # nonlinear_problem = NonlinearVariationalProblem(self.problem.F, self.problem.up_next, self.problem.bd.bcs, J)
        # nonlinear_solver  = NonlinearVariationalSolver(nonlinear_problem)

        # ### Set some parameters ###
        # solver_parameters = nonlinear_solver.parameters
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


        nonlinear_solver = self.params["solver"].get("nonlinear_solver", "snes")
        relaxation = self.params["solver"].get("newton_relaxation", 1.0)

        self.fprint("Solving with {0}".format(nonlinear_solver))
        if nonlinear_solver == "newton":
            self.fprint("Relaxation parameter = {: 1.2f}".format(relaxation))

            newton_options = {"relaxation_parameter": relaxation,
                              "maximum_iterations": 40,
                              "linear_solver": "mumps",
                              "absolute_tolerance": 1e-6,
                              "relative_tolerance": 1e-5}
        
            solver_parameters = {"nonlinear_solver": "newton",
                                 "newton_solver": newton_options}

        elif nonlinear_solver == "snes":
            # ### Add some helper functions to solver options ###
            solver_parameters = {"nonlinear_solver": "snes",
                                 "snes_solver": {
                                 "linear_solver": "mumps", 
                                 "maximum_iterations": 40,
                                 "error_on_nonconvergence": True,
                                 "line_search": "bt",
                                 }}
                                 
        else:
            raise ValueError("Unknown nonlinear solver type: {0}".format(nonlinear_solver))

        ### Start the Solve Process ###
        self.fprint("Solving",special="header")
        start = time.time()
        
        # ### Solve the Baseline Problem ###
        # solve(self.problem.F_sans_tf == 0, self.problem.up_next, self.problem.bd.bcs, solver_parameters=solver_parameters, **self.extra_kwarg)

        # ### Store the Baseline and Assign for the real solve ###
        # self.up_baseline = self.problem.up_next.copy(deepcopy=True)
        # self.problem.up_next.assign(self.up_baseline)

        ### Solve the real problem ###
        solve(self.problem.F == 0, self.problem.up_next, self.problem.bd.bcs, solver_parameters=solver_parameters)
        stop = time.time()
        self.fprint("Solve Complete: {:1.2f} s".format(stop-start),special="footer")
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
        if self.optimizing or self.save_power:
            self.J += -self.CalculatePowerFunctional((iter_val-self.problem.dom.init_wind)) 

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
        self.angles = np.linspace(self.wind_range[0],self.wind_range[1],self.num_wind,endpoint=self.endpoint)

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