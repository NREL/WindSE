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

        self.control_types = self.params["optimization"]["controls"]
        self.layout_bounds = self.params["optimization"].get("layout_bounds",[self.farm.ex_x,self.farm.ex_y])

        self.iteration = 0

        self.fprint("Setting Up Optimizer",special="header")
        
        self.fprint("Controls: {0}".format(self.control_types))
        self.CreateControls()
 
        
        self.fprint("Define Bounds")
        self.CreateBounds()

        self.get_minimum_distance_constraint_func(self.controls, 200)

        self.fprint("Define Optimizing Functional")
        if self.params["solver"]["type"] == "multiangle":
            self.J = self.solver.J
            self.Jhat = ReducedFunctional(self.J, self.controls, eval_cb_post=self.ReducedFunctionalCallback)
            self.Jcurrent = self.J
        else:
            self.PowerFunctional()

        self.fprint("Number of Controls: {:d}".format(len(self.controls)),special="header")
        self.OptPrintFunction(self.init_vals)
        self.fprint("",special="footer")


        self.fprint("Optimizer Setup",special="footer")

    def ReducedFunctionalCallback(self, j, m):
        self.Jcurrent = j 

    def CreateControls(self):
        self.controls = []
        self.names = []
        self.indexes = [[],[],[],[]]
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

    def CreateBounds(self):
        lower_bounds = []
        upper_bounds = []

        if "layout" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant((self.layout_bounds[0][0] + self.farm.radius[i])))
                lower_bounds.append(Constant((self.layout_bounds[1][0] + self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[0][1] - self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[1][1] - self.farm.radius[i])))

        if "yaw" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(-pi/4.))
                upper_bounds.append(Constant(pi/4.))

        if "axial" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(1.))

        self.bounds = [lower_bounds,upper_bounds]

    def AssignControls(self):
        if "layout" in self.control_types:
            self.farm.mx = self.controls[self.indexes[0]]
            self.farm.my = self.controls[self.indexes[1]]
        if "yaw" in self.control_types:
            self.farm.myaw = self.controls[self.indexes[2]]
        if "axial" in self.control_types:
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

        if self.farm.yaw[0]**2 > 1e-4:
            self.J = assemble(-dot(self.problem.tf,self.solver.u_next),*dx)
        else:
            self.J = assemble(-inner(dot(self.problem.tf,self.solver.u_next),self.solver.u_next[0]**2+self.solver.u_next[1]**2)*dx)
        self.Jhat = ReducedFunctional(self.J, self.controls, eval_cb_post=self.ReducedFunctionalCallback) 
        self.Jcurrent = self.J

    def ListControls(self,m):
        for i,val in enumerate(m):
            self.fprint(self.names[i]+": "+repr(float(m[i])))

    def PlotLayout(self,m,show=False):
        x_val = []
        y_val = []

        for i,val in enumerate(m):
            if "x" in self.names[i]:
                x_val.append(float(m[i]))
            elif "y" in self.names[i]:
                y_val.append(float(m[i]))

        z_val = self.problem.dom.Ground(x_val,y_val)

        ### Create the path names ###
        folder_string = self.params.folder+"/plots/"
        file_string = self.params.folder+"/plots/wind_farm_step_"+repr(self.iteration)+".pdf"

        ### Check if folder exists ###
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        ### Create a list that outlines the extent of the farm ###
        ex_list_x = [self.layout_bounds[0][0],self.layout_bounds[0][1],self.layout_bounds[0][1],self.layout_bounds[0][0],self.layout_bounds[0][0]]
        ex_list_y = [self.layout_bounds[1][0],self.layout_bounds[1][0],self.layout_bounds[1][1],self.layout_bounds[1][1],self.layout_bounds[1][0]]

        ### Generate and Save Plot ###
        plt.clf()
        if hasattr(self.problem.dom,"boundary_line"):
            plt.plot(*self.problem.dom.boundary_line,c="k")
        plt.plot(ex_list_x,ex_list_y,c="r")
        p=plt.scatter(x_val,y_val,c=range(self.farm.numturbs))
        plt.xlim(self.problem.dom.x_range[0],self.problem.dom.x_range[1])
        plt.ylim(self.problem.dom.y_range[0],self.problem.dom.y_range[1])
        clb = plt.colorbar(p)
        clb.ax.set_ylabel('Hub Height')

        plt.title("Power Output: {: 5.2f}".format(self.Jcurrent))
        plt.savefig(file_string)
        if show:
            plt.show()

    def OptPrintFunction(self,m):
        self.ListControls(m)

        if "layout" in self.control_types:
            self.PlotLayout(m,show=False)
        
        self.iteration += 1

    def get_minimum_distance_constraint_func(self, m_pos, min_distance=200):
        if "layout" in self.control_types and len(self.control_types)==1:
            self.dist_constraint = MinimumDistanceConstraint(m_pos, min_distance)
        else:
            print("minimum distance is supported when only optimizing layout")
            self.dist_constraint = None

    def Optimize(self):

        self.fprint("Beginning Optimization",special="header")

        m_opt=minimize(self.Jhat, method="SLSQP", options = {"disp": True}, constraints = self.dist_constraint, bounds = self.bounds, callback = self.OptPrintFunction)
        # m_opt=minimize(self.Jhat, method="L-BFGS-B", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)

        self.fprint("Assigning New Values")
        self.AssignControls()
        
        self.fprint("Solving With New Values")
        self.solver.Solve()

        self.fprint("Optimization Finished",special="header")

        return m_opt

    def TaylorTest(self):
        
        self.fprint("Beginning Taylor Test",special="header")

        # h = [Constant(10)]*(len(self.controls))
        h = [Constant(0.001)]*(len(self.controls))

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
        print("In mimimum distance constraint")

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

        print("In mimimum distance constraint function eval")
        # print "distances: ", arr*lengthscale
        numClose = 0
        for i in range(len(arr)):
            if arr[i]<0:
                print(arr[i]*lengthscale)
                numClose +=1
        print("Number of turbines in violation of spacing constraint:", numClose)
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
        print("In mimimum distance constraint Jacobian eval")
        return np.array(ieqcons)