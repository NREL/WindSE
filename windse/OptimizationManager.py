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
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np
    import copy
    from sys import platform
    import time
    from memory_profiler import memory_usage

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
else:
    InequalityConstraint = object

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
        self.Sx = self.problem.dom.xscale

        self.control_types = self.params["optimization"]["controls"]
        self.layout_bounds = np.array(self.params["optimization"].get("layout_bounds",[np.array(self.farm.ex_x)/self.Sx,np.array(self.farm.ex_y)/self.Sx]))*self.Sx

        self.iteration = 0

        self.fprint("Setting Up Optimizer",special="header")
        
        self.fprint("Controls: {0}".format(self.control_types))
        self.CreateControls()
 
        
        self.fprint("Define Bounds")
        self.CreateBounds()

        self.get_minimum_distance_constraint_func(self.controls, 2*np.mean(self.problem.farm.HH))

        self.fprint("Define Optimizing Functional")
        self.J = self.solver.J
        self.Jhat = ReducedFunctional(self.J, self.controls, eval_cb_post=self.ReducedFunctionalCallback)
        self.Jcurrent = self.J

        self.fprint("Number of Controls: {:d}".format(len(self.controls)),special="header")
        self.OptPrintFunction(self.init_vals)
        self.fprint("",special="footer")


        self.fprint("Optimizer Setup",special="footer")

    def ReducedFunctionalCallback(self, j, m):
        self.Jcurrent = j 

    def CreateControls(self):
        self.controls = []
        self.names = []
        self.indexes = [[],[],[],[],[],[]]
        self.init_vals = []
        j = 0
        if "layout" in self.control_types:
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

        if "yaw" in self.control_types:
            for i in range(self.farm.numturbs):
                self.indexes[2].append(j)
                j+=1
                self.names.append("yaw_"+repr(i))
                self.controls.append(Control(self.farm.myaw[i]))
                self.init_vals.append(self.farm.myaw[i])

        if "axial" in self.control_types:
            for i in range(self.farm.numturbs):
                self.indexes[3].append(j)
                j+=1
                self.names.append("axial_"+repr(i))
                self.controls.append(Control(self.farm.ma[i]))
                self.init_vals.append(self.farm.ma[i])

        if "lift" in self.control_types:
            for i in range(self.problem.num_blade_segments):
                self.indexes[4].append(j)
                j+=1
                self.names.append("lift_"+repr(i))
                self.controls.append(Control(self.problem.mcl[i]))
                self.init_vals.append(self.problem.mcl[i])

        if "drag" in self.control_types:
            for i in range(self.problem.num_blade_segments):
                self.indexes[5].append(j)
                j+=1
                self.names.append("drag_"+repr(i))
                self.controls.append(Control(self.problem.mcd[i]))
                self.init_vals.append(self.problem.mcl[i])

    def CreateBounds(self):
        lower_bounds = []
        upper_bounds = []

        if "layout" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant((self.layout_bounds[0][0])))# + self.farm.radius[i])))
                lower_bounds.append(Constant((self.layout_bounds[1][0])))# + self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[0][1])))# - self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[1][1])))# - self.farm.radius[i])))

        if "yaw" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(0.0))
                upper_bounds.append(Constant(25*pi/180.0))

        if "axial" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(1.))

        if "lift" in self.control_types:
            for i in range(self.problem.num_blade_segments):
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(2.))

        if "drag" in self.control_types:
            for i in range(self.problem.num_blade_segments):
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(2.))

        self.bounds = [lower_bounds,upper_bounds]

    def ListControls(self,m):

        self.fprint("Current Objective Value: " + repr(float(self.Jcurrent)))

        if "layout" in self.control_types:
            for i in range(self.farm.numturbs):
                self.fprint("Location Turbine {0:} of {1:}: {2: 4.2f}, {3: 4.2f}".format(i+1,self.farm.numturbs,self.farm.x[i],self.farm.y[i]))

        if "yaw" in self.control_types:
            for i in range(self.farm.numturbs):
                self.fprint("Yaw Turbine {0:} of {1:}: {2: 4.6f}".format(i+1,self.farm.numturbs,self.farm.yaw[i]))

        for i, val in enumerate(m):
            self.fprint(self.names[i] +": " +repr(float(val)))

    def SaveControls(self,m):

        folder_string = self.params.folder+"/data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)


        new_values = {}
        m_f = np.array(m,dtype=float)
        if "layout" in self.control_types:
            new_values["x"]   = m_f[self.indexes[0]]
            new_values["y"]   = m_f[self.indexes[1]]
        if "yaw" in self.control_types:
            new_values["yaw"] = m_f[self.indexes[2]]
        if "axial" in self.control_types:
            new_values["a"]   = m_f[self.indexes[3]]

        self.problem.farm.UpdateControls(**new_values)
        self.problem.farm.SaveWindFarm(val=self.iteration)

        if self.iteration == 0:

            #### ADD HEADER ####
            self.last_m = np.zeros(len(m))
            for i in range(len(m)):
                self.last_m[i]=float(m[i])
            err = 0.0
            f = open(folder_string+"opt_data.txt",'wb')
        else:
            err = np.linalg.norm(m-self.last_m)
            self.last_m = m
            f = open(folder_string+"opt_data.txt",'ab')

        output_data = np.concatenate(((self.Jcurrent, err),m))

        np.savetxt(f,[output_data])
        f.close()

    def OptPrintFunction(self,m):
        # tape = get_working_tape()
        # tape.visualise()
        # exit()



        mem0=memory_usage()[0]
        tick = time.time()
        mem_out, der = memory_usage(self.Jhat.derivative,max_usage=True,retval=True,max_iterations=1)
        for i, d in enumerate(der):
            self.fprint("dJd"+self.names[i] +": " +repr(float(d)))
        tock = time.time()
        self.fprint("Time Elapsed: {:1.2f} s".format(tock-tick))
        self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))
        # exit()

        # exit()
        # if "layout" in self.control_types:
        #     self.PlotLayout(m,show=False)

        self.SaveControls(m)
        self.ListControls(m)

        self.problem.farm.Plot(filename="wind_farm_step_"+repr(self.iteration),power=-self.Jcurrent)
        
        self.iteration += 1

    def get_minimum_distance_constraint_func(self, m_pos, min_distance=200):
        if "layout" in self.control_types and len(self.control_types)==1:
            self.dist_constraint = MinimumDistanceConstraint(m_pos, min_distance)
        else:
            print("minimum distance is supported when only optimizing layout")
            self.dist_constraint = None

    def Optimize(self):

        self.fprint("Beginning Optimization",special="header")




        if "layout" in self.control_types:
            m_opt=minimize(self.Jhat, method="SLSQP", options = {"disp": True}, constraints = self.dist_constraint, bounds = self.bounds, callback = self.OptPrintFunction)
        else:
            m_opt=minimize(self.Jhat, method="SLSQP", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)
        # m_opt=minimize(self.Jhat, method="L-BFGS-B", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)

        self.fprint("Assigning New Values")
        new_values = {}
        m_f = np.array(m_opt,dtype=float)
        if "layout" in self.control_types:
            new_values["x"]   = m_f[self.indexes[0]]
            new_values["y"]   = m_f[self.indexes[1]]
        if "yaw" in self.control_types:
            new_values["yaw"] = m_f[self.indexes[2]]
        if "axial" in self.control_types:
            new_values["a"]   = m_f[self.indexes[3]]
        self.problem.farm.UpdateControls(**new_values)
        # self.AssignControls()
        
        self.fprint("Solving With New Values")
        self.solver.Solve()

        self.fprint("Optimization Finished",special="header")

        return m_opt

    def TaylorTest(self):
        
        self.fprint("Beginning Taylor Test",special="header")

        h = []
        for i,c in enumerate(self.controls):
            h.append(Constant(1*max(abs(float(self.bounds[1][i])),abs(float(self.bounds[1][i])))))

        conv_rate = taylor_test(self.Jhat, self.init_vals, h)

        self.fprint("Convergence Rates:")
        self.fprint("")
        self.fprint(conv_rate)
        self.fprint("")

        self.fprint("Taylor Test Finished",special="footer")



        return conv_rate

class MinimumDistanceConstraint(InequalityConstraint):
    def __init__(self, m_pos, min_distance=200):

        self.min_distance = min_distance
        self.m_pos = m_pos
        # print("In mimimum distance constraint")

    def length(self):
        nconstraints = comb(len(self.m_pos)/2,2.)
        return nconstraints

    def function(self, m):
        ieqcons = []
        
        m_pos = m

        for i in range(int(len(m_pos) / 2)):
            for j in range(int(len(m_pos) / 2)):
                if j > i:
                    ieqcons.append(((m_pos[2 * i] - m_pos[2 * j])**2 + (m_pos[2 * i + 1] - m_pos[2 * j + 1])**2) - self.min_distance**2)

        arr = np.array(ieqcons)

        # print("In mimimum distance constraint function eval")
        # print "distances: ", arr*lengthscale
        numClose = 0
        for i in range(len(arr)):
            if arr[i]<0:
                # print(arr[i]*lengthscale)
                numClose +=1
        if numClose > 1:
            print("Warning: Number of turbines in violation of spacing constraint: "+repr(numClose))
        return np.array(ieqcons)

    def jacobian(self, m):
        ieqcons = []
        
        m_pos = m

        for i in range(int(len(m_pos) / 2)):
            for j in range(int(len(m_pos) / 2)):
                if j>i:
                    prime_ieqcons = np.zeros(len(m))

                    prime_ieqcons[2 * i] = 2 * (m_pos[2 * i] - m_pos[2 * j])
                    prime_ieqcons[2 * j] = -2 * (m_pos[2 * i] - m_pos[2 * j])
                    prime_ieqcons[2 * i + 1] = 2 * (m_pos[2 * i + 1] - m_pos[2 * j + 1])
                    prime_ieqcons[2 * j + 1] = -2 * (m_pos[2 * i + 1] - m_pos[2 * j + 1])

                    ieqcons.append(prime_ieqcons)
        # print("In mimimum distance constraint Jacobian eval")
        return np.array(ieqcons)