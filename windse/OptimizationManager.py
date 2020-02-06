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
import numpy as np

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

        self.control_types = self.params["optimization"]["controls"]
        print('control type:', self.control_types)
        self.layout_bounds = self.params["optimization"].get("layout_bounds",[self.farm.ex_x,self.farm.ex_y])

        self.iteration = 0

        self.fprint("Setting Up Optimizer",special="header")

        self.fprint("Controls: {0}".format(self.control_types))
        self.CreateControls()


        self.fprint("Define Bounds")
        self.CreateBounds()

        self.get_minimum_distance_constraint_func(self.controls, 2*np.mean(self.problem.farm.HH))

        self.fprint("Define Optimizing Functional")
        if self.params["solver"]["type"] == "multiangle":
            self.J = assemble(self.solver.J)
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
                lower_bounds.append(Constant((self.layout_bounds[0][0])))# + self.farm.radius[i])))
                lower_bounds.append(Constant((self.layout_bounds[1][0])))# + self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[0][1])))# - self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[1][1])))# - self.farm.radius[i])))

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

        self.J = assemble(-(dot(self.problem.tf,self.solver.u_next))*dx)
        self.Jhat = ReducedFunctional(self.J, self.controls, eval_cb_post=self.ReducedFunctionalCallback)
        self.Jcurrent = self.J

    def Gradient(self):
        """
        Returns a gradient of the objective function
        """
        dJdma= compute_gradient(self.J, self.controls)
        # gradient_list = [float(dd) for dd in dJdma] # Gradient given in a list from Fenics
        gradient_list = np.array(dJdma, dtype=np.float)
        print(gradient_list)

        return gradient_list

    def ListControls(self,m):
        if "layout" in self.control_types:
            for i in range(self.farm.numturbs):
                self.fprint("Turbine {0:} of {1:}: {2: 4.2f}, {3: 4.2f}".format(i,self.farm.numturbs,self.x_val[i],self.y_val[i]))

    def PlotLayout(self,m,show=False):
        self.x_val = []
        self.y_val = []

        for i,val in enumerate(m):
            if "x_" in self.names[i]:
                self.x_val.append(float(m[i]))
            elif "y_" in self.names[i]:
                self.y_val.append(float(m[i]))

        z_val = self.problem.dom.Ground(self.x_val,self.y_val)+self.problem.farm.HH

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
        p=plt.scatter(self.x_val,self.y_val,c=z_val,cmap="coolwarm",edgecolors=(0, 0, 0, 1))
        # p=plt.scatter(self.x_val,self.y_val,c="k",s=70)
        # p=plt.scatter(self.x_val,self.y_val,c=range(self.farm.numturbs))
        plt.xlim(self.problem.dom.x_range[0],self.problem.dom.x_range[1])
        plt.ylim(self.problem.dom.y_range[0],self.problem.dom.y_range[1])
        clb = plt.colorbar(p)
        clb.ax.set_ylabel('Hub Elevation')

        plt.title("Power Output: {: 5.2f}".format(-self.Jcurrent))
        plt.savefig(file_string, transparent=True)
        if show:
            plt.show()

    def SaveControls(self,m):

        folder_string = self.params.folder+"/data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        if "layout" in self.control_types:
            self.problem.farm.UpdateConstants(m=m,control_types=self.control_types,indexes=self.indexes)
            self.problem.farm.SaveWindFarm(val=self.iteration)

        if self.iteration == 0:
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

        if "layout" in self.control_types:
            self.PlotLayout(m,show=False)

        self.SaveControls(m)
        self.ListControls(m)

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
            # m_opt=minimize(self.Jhat, method="SLSQP", options = {"disp": True}, constraints = self.dist_constraint, bounds = self.bounds, callback = self.OptPrintFunction)
            m_opt=minimize(self.Jhat, method="L-BFGS-B", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)
        else:
            m_opt=minimize(self.Jhat, method="SLSQP", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)
        # m_opt=minimize(self.Jhat, method="L-BFGS-B", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)

        self.fprint("Assigning New Values")
        self.problem.farm.UpdateConstants(m=m_opt,control_types=self.control_types,indexes=self.indexes)
        # self.AssignControls()

        self.fprint("Solving With New Values")
        self.solver.Solve()

        self.fprint("Optimization Finished",special="header")

        return m_opt

    def TaylorTest(self):

        self.fprint("Beginning Taylor Test",special="header")

        h = [Constant(10)]*(len(self.controls))
        # h = [Constant(0.001)]*(len(self.controls))

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
