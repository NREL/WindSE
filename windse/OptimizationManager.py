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
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"
    
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
    if windse_parameters.dolfin_adjoint:
        from dolfin_adjoint import *

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt
else:
    InequalityConstraint = object
    
    
import openmdao.api as om


class ObjComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('initial_DVs', types=np.ndarray)
        self.options.declare('J', types=object)
        self.options.declare('dJ', types=object)
        self.options.declare('callback', types=object)
        
    def setup(self):
        self.add_input('DVs', val=self.options['initial_DVs'])
        self.add_output('obj', val=0.)

        self.declare_partials('*', '*')
        
    def compute(self, inputs, outputs):
        m = list(inputs['DVs'])
        computed_output = self.options['J'](m)
        outputs['obj'] = computed_output
        self.options['callback'](m)
        
    def compute_partials(self, inputs, partials):
        m = list(inputs['DVs'])
        jac = self.options['dJ'](m)
        partials['obj', 'DVs'] = jac
        
        
class ConsComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('initial_DVs', types=np.ndarray)
        self.options.declare('J', types=object)
        self.options.declare('dJ', types=object)
        self.options.declare('con_name', types=str)
        
    def setup(self):
        self.con_name = self.options["con_name"]
        self.add_input('DVs', val=self.options['initial_DVs'])
        
        output = self.options['J'](self.options['initial_DVs'])
        self.add_output(self.con_name, val=output)

        self.declare_partials('*', '*')
        
    def compute(self, inputs, outputs):
        m = list(inputs['DVs'])
        computed_output = self.options['J'](m)
        outputs[self.con_name] = computed_output
        
    def compute_partials(self, inputs, partials):
        m = list(inputs['DVs'])
        jac = self.options['dJ'](m)
        partials[self.con_name, 'DVs'] = jac
        

def gather(m):
    if isinstance(m, list):
        return list(map(gather, m))
    elif hasattr(m, "_ad_to_list"):
        return m._ad_to_list(m)
    else:
        return m  # Assume it is gathered already
    
def om_wrapper(J, initial_DVs, dJ, H, bounds, **kwargs):
    
    # build the model
    prob = om.Problem(model=om.Group())
    
    if 'callback' in kwargs:
        callback = kwargs['callback']
    else:
        callback = None

    prob.model.add_subsystem('obj_comp', ObjComp(initial_DVs=initial_DVs, J=J, dJ=dJ, callback=callback), promotes=['*'])
    
    constraint_types = []
    if 'constraints' in kwargs:
        constraints = kwargs['constraints']
        
        if not isinstance(constraints, list):
            constraints = [constraints]
        
        for idx, c in enumerate(constraints):
            if isinstance(c, InequalityConstraint):
                typestr = "ineq"
            elif isinstance(c, EqualityConstraint):
                typestr = "eq"
            else:
                raise Exception("Unknown constraint class")

            def jac(x):
                out = c.jacobian(x)
                return [gather(y) for y in out]

            constraint_types.append(typestr)
            
            con_name = f'con_{idx}'
            prob.model.add_subsystem(f'cons_comp_{idx}', ConsComp(initial_DVs=initial_DVs, J=c.function, dJ=jac, con_name=con_name), promotes=['*'])
    
    lower_bounds = bounds[:, 0]
    upper_bounds = bounds[:, 1]
    
    # set up the optimization
    if 'SLSQP' in kwargs['opt_routine']:
        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['optimizer'] = 'SLSQP'
    elif 'SNOPT' in kwargs['opt_routine']:
        prob.driver = om.pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SNOPT'
    
    prob.model.add_design_var('DVs', lower=lower_bounds, upper=upper_bounds)
    prob.model.add_objective('obj')
    
    for idx, constraint_type in enumerate(constraint_types):
        con_name = f'con_{idx}'
        if constraint_type == "eq":
            prob.model.add_constraint(con_name, equals=0.)
        else:
            # Inequality means it's positive from scipy and dolfin
            prob.model.add_constraint(con_name, lower=0.)            

    prob.setup()
    
    prob.set_val('DVs', initial_DVs)

    # Run the optimization
    prob.run_driver()
    
    # Return the optimal design variables
    return(prob['DVs'])
    

    
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
        self.tag_output = self.params.tag_output
        self.debug_mode = self.params.debug_mode
        self.xscale = self.problem.dom.xscale

        ### Update attributes based on params file ###
        for key, value in self.params["optimization"].items():
            if isinstance(value,list):
                setattr(self,key,np.array(value))
            else:
                setattr(self,key,value)

        ### Process parameters ###
        if "layout" in self.control_types:
            if isinstance(self.layout_bounds,(list, np.ndarray)):
                self.layout_bounds = np.array(self.layout_bounds)
            elif self.layout_bounds == "wind_farm":
                self.layout_bounds = np.array([self.farm.ex_x,self.farm.ex_y])
            else:
                self.layout_bounds = np.array([[0,0],[0,0]])
            self.layout_bounds = self.layout_bounds*self.xscale

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
        self.OptPrintFunction(self.init_vals,None)
        self.fprint("",special="footer")

        self.fprint("Optimizer Setup",special="footer")

    def DebugOutput(self):
        if self.debug_mode:

            self.tag_output("n_controls", len(self.controls))
            self.tag_output("obj_value", float(self.J))

            ### Output initial control values ###
            for i, val in enumerate(self.controls):
                self.tag_output("val0_"+self.names[i],val.values())

            ### Output gradient ###
            if hasattr(self,"gradient"):
                for i, d in enumerate(self.gradients):
                    self.tag_output("grad_"+self.names[i],float(d))
            
            ### TODO: Output taylor convergence data
            if hasattr(self,"conv_rate"):
                pass

            ### TODO: Output optimized controls
            if hasattr(self,"m_opt"):
                pass


    def RecomputeReducedFunctional(self):
        self.CreateControls()
        self.J = self.solver.J
        self.Jhat = ReducedFunctional(self.J, self.controls, eval_cb_post=self.ReducedFunctionalCallback)

        self.Jcurrent = self.J

    def ReducedFunctionalCallback(self, j, m):
        self.Jcurrent = j 

    def CreateControls(self):

        ### Control pointers is a way of reference which parts of the original farm lists correspond to the controls. it is rather inefficient but it works so...
        self.controls = []
        self.control_pointers = []
        self.names = []
        self.indexes = [[],[],[],[],[],[],[]]
        self.init_vals = []
        j = 0
        if "layout" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[0].append(j)
                j+=1
                self.names.append("x_"+repr(i))
                self.controls.append(Control(self.farm.mx[i]))
                self.control_pointers.append((self.farm.x,i))
                self.init_vals.append(self.farm.mx[i])

                self.indexes[1].append(j)
                j+=1
                self.names.append("y_"+repr(i))
                self.controls.append(Control(self.farm.my[i]))
                self.control_pointers.append((self.farm.y,i))
                self.init_vals.append(self.farm.my[i])

        if "yaw" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[2].append(j)
                j+=1
                self.names.append("yaw_"+repr(i))
                self.controls.append(Control(self.farm.myaw[i]))
                self.control_pointers.append((self.farm.yaw,i))
                self.init_vals.append(self.farm.myaw[i])

        if "axial" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[3].append(j)
                j+=1
                self.names.append("axial_"+repr(i))
                self.controls.append(Control(self.farm.ma[i]))
                self.control_pointers.append((self.farm.a,i))
                self.init_vals.append(self.farm.ma[i])

        if "lift" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    self.control_pointers.append((self.farm.cl,[i,k]))
                    self.indexes[4].append(j)
                    j+=1
                    self.names.append("lift_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.mcl[i][k]))
                    self.init_vals.append(self.farm.mcl[i][k])

        if "drag" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    self.control_pointers.append((self.farm.cd,[i,k]))
                    self.indexes[5].append(j)
                    j+=1
                    self.names.append("drag_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.mcd[i][k]))
                    self.init_vals.append(self.farm.mcd[i][k])

        if "chord" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    self.control_pointers.append((self.farm.chord,[i,k]))
                    self.indexes[6].append(j)
                    j+=1
                    self.names.append("chord_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.mchord[i][k]))
                    self.init_vals.append(self.farm.mchord[i][k])
        self.num_controls = len(self.controls)

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
                lower_bounds.append(Constant(-45*pi/180.0))
                upper_bounds.append(Constant(45*pi/180.0))

        if "axial" in self.control_types:
            for i in range(self.farm.numturbs):
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(1.))

        if "lift" in self.control_types:
            for i in self.solver.opt_turb_id:
                for i in range(self.farm.num_blade_segments):
                    lower_bounds.append(Constant(0))
                    upper_bounds.append(Constant(2.))

        if "drag" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    lower_bounds.append(Constant(0))
                    upper_bounds.append(Constant(2.))

        if "chord" in self.control_types:
            for i in self.solver.opt_turb_id:
                c_avg = 0
                for k in range(self.farm.num_blade_segments):
                    modifier = 2.0
                    max_chord = self.farm.max_chord
                    seg_chord = self.farm.baseline_chord[k]
                    lower_bounds.append(Constant(seg_chord/modifier))
                    upper_bounds.append(Constant(np.maximum(np.minimum(seg_chord*modifier,max_chord),c_avg)))
                    c_avg = (c_avg*k+seg_chord)/(k+1)

        self.bounds = [lower_bounds,upper_bounds]

    def Gradient(self):
        """
        Returns a gradient of the objective function
        """
        mem0=memory_usage()[0]
        tick = time.time()

        capture_memory = False
        if capture_memory:
            mem_out, der = memory_usage(self.Jhat.derivative,max_usage=True,retval=True,max_iterations=1)
        else: 
            mem_out = 2*mem0
            der = self.Jhat.derivative()

        folder_string = self.params.folder+"data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)
        f = open(folder_string+"gradient_data.txt",'w')
        f_header = "control    value    derivative"
        f.write(f_header+"\n")

        for i, d in enumerate(der):
            ctl_val = float(self.controls[i].values())
            d_out = str(self.names[i] + "    " +repr(ctl_val)+ "    " +repr(float(d)))
            self.fprint(d_out)
            f.write(d_out+"\n")
        f.close()

        tock = time.time()
        self.fprint("Time Elapsed: {:1.2f} s".format(tock-tick))
        if capture_memory:
            self.fprint("Memory Used:  {:1.2f} MB".format(mem_out-mem0))

        self.gradients = np.array(der, dtype=float)
        self.DebugOutput()

        return self.gradients

    def ListControls(self,m):
        self.fprint("Iteration "+repr(self.iteration)+" Complete",special="header")
        self.fprint("Current Objective Value: " + repr(float(self.Jcurrent)))

        # if "layout" in self.control_types:
        #     for i in range(self.farm.numturbs):
        #         self.fprint("Location Turbine {0:} of {1:}: {2: 4.2f}, {3: 4.2f}".format(i+1,self.farm.numturbs,self.farm.x[i],self.farm.y[i]))

        # if "yaw" in self.control_types:
        #     for i in range(self.farm.numturbs):
        #         self.fprint("Yaw Turbine {0:} of {1:}: {2: 4.6f}".format(i+1,self.farm.numturbs,self.farm.yaw[i]))
        self.fprint("Previous Control Values",special="header")
        for i, [l, ix] in enumerate(self.control_pointers):
            if not isinstance(ix,int):
                self.fprint(self.names[i] +": " +repr(float(l[ix[0]][ix[1]])))
            else:
                self.fprint(self.names[i] +": " +repr(float(l[ix])))
        self.fprint("",special="footer")

        self.fprint("Next Control Values",special="header")
        for i, val in enumerate(m):
            self.fprint(self.names[i] +": " +repr(float(val)))
        self.fprint("",special="footer")

        self.fprint("Iteration "+repr(self.iteration)+" Complete",special="footer")

    def SaveControls(self,m):

        folder_string = self.params.folder+"/data/"
        if not os.path.exists(folder_string): os.makedirs(folder_string)

        # new_values = {}
        # m_f = np.array(self.control_pointers,dtype=float)
        # if "layout" in self.control_types:
        #     new_values["x"]   = m_f[self.indexes[0]]
        #     new_values["y"]   = m_f[self.indexes[1]]
        # if "yaw" in self.control_types:
        #     new_values["yaw"] = m_f[self.indexes[2]]
        # if "axial" in self.control_types:
        #     new_values["a"]   = m_f[self.indexes[3]]

        # self.problem.farm.UpdateControls(**new_values)
        self.problem.farm.SaveWindFarm(val=self.iteration)
        # print(m)
        # print(type(m))
        m_new = np.array(m,dtype=float)
        m_old = []
        for l, i in self.control_pointers:
            if not isinstance(i,int):
                m_old.append(float(l[i[0]][i[1]]))
            else:
                m_old.append(float(l[i]))
        # print(m_new)
        # print(type(m_new))            
        if self.iteration == 0:

            #### ADD HEADER ####
            self.last_m = np.zeros(self.num_controls)
            for i in range(self.num_controls):
                self.last_m[i]=float(m_new[i])
            err = 0.0
            f = open(folder_string+"opt_data.txt",'w')
            header = str("Objective    Change    Prev_Controls:    p_"+"    p_".join(self.names)+"    New_Controls:    n_"+"    n_".join(self.names)+"\n")
            f.write(header)
        else:
            err = np.linalg.norm(m_new-self.last_m)
            self.last_m = copy.copy(m_new)
            f = open(folder_string+"opt_data.txt",'a')

        output_data = np.concatenate(((self.Jcurrent, err, self.num_controls),m_old))
        output_data = np.concatenate((output_data,(self.num_controls,)))
        output_data = np.concatenate((output_data,m_new))

        np.savetxt(f,[output_data])
        f.close()

    def OptPrintFunction(self,m,test=None):
        if test is not None:
            print("Hey, this method actually gives us more info")
        # print(np.array(m,dtype=float))
        # print(np.array(self.control_pointers,dtype=float))
        # print(np.array(self.problem.farm.myaw,dtype=float))
        self.SaveControls(m)
        self.ListControls(m)

        if "layout" in self.control_types or "yaw" in self.control_types:
            self.problem.farm.PlotFarm(filename="wind_farm_step_"+repr(self.iteration),power=-self.Jcurrent)

        if "chord" in self.control_types:
            c_lower = np.array(self.bounds[0])[self.indexes[6]] 
            c_upper = np.array(self.bounds[1])[self.indexes[6]] 
            self.problem.farm.PlotChord(filename="chord_step_"+repr(self.iteration),power=-self.Jcurrent,bounds=[c_lower,c_upper])
        
        self.iteration += 1

    def get_minimum_distance_constraint_func(self, m_pos, min_distance=200):
        if "layout" in self.control_types and len(self.control_types)==1:
            self.dist_constraint = MinimumDistanceConstraint(m_pos, min_distance)
        else:
            print("minimum distance is supported when only optimizing layout")
            self.dist_constraint = None

    def Optimize(self):

        self.fprint("Beginning Optimization",special="header")
        
        # TODO : simplify this logic
        if "SNOPT" in self.opt_routine or "OM_SLSQP" in self.opt_routine:
            if "layout" in self.control_types:
                m_opt=minimize(self.Jhat, method="Custom", options = {"disp": True}, constraints = self.dist_constraint, bounds = self.bounds, callback = self.OptPrintFunction, algorithm=om_wrapper, opt_routine=self.opt_routine)
            else:
                m_opt=minimize(self.Jhat, method="Custom", options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction, algorithm=om_wrapper, opt_routine=self.opt_routine)
        else:
            if "layout" in self.control_types:
                m_opt=minimize(self.Jhat, method=self.opt_routine, options = {"disp": True}, constraints = self.dist_constraint, bounds = self.bounds, callback = self.OptPrintFunction)
            else:
                m_opt=minimize(self.Jhat, method=self.opt_routine, options = {"disp": True}, bounds = self.bounds, callback = self.OptPrintFunction)
                
        self.m_opt = m_opt

        if self.num_controls == 1:
            self.m_opt = (self.m_opt,)
        self.OptPrintFunction(m_opt)
        # self.fprint("Assigning New Values")
        # new_values = {}
        # m_f = np.array(m_opt,dtype=float)
        # if "layout" in self.control_types:
        #     new_values["x"]   = m_f[self.indexes[0]]
        #     new_values["y"]   = m_f[self.indexes[1]]
        # if "yaw" in self.control_types:
        #     new_values["yaw"] = m_f[self.indexes[2]]
        # if "axial" in self.control_types:
        #     new_values["a"]   = m_f[self.indexes[3]]
        # self.problem.farm.UpdateControls(**new_values)
        
        # self.fprint("Solving With New Values")
        # self.solver.Solve()
        self.DebugOutput()
        self.fprint("Optimization Finished",special="footer")

        return self.m_opt

    def TaylorTest(self):
        
        self.fprint("Beginning Taylor Test",special="header")

        h = []
        for i,c in enumerate(self.controls):
            h.append(Constant(0.1))
            # h.append(Constant(0.01*max(abs(float(self.bounds[1][i])),abs(float(self.bounds[1][i])))))
            # h.append(Constant(10.0*abs(float(self.bounds[1][i])-float(self.bounds[0][i]))/2.0))
            # h.append(Constant(0.01*abs(np.mean(self.bounds[1])+np.mean(self.bounds[0]))/2.0))

        print(np.array(h,dtype=float))

        self.conv_rate = taylor_test(self.Jhat, self.init_vals, h)
        self.DebugOutput()

        self.fprint("Convergence Rates:")
        self.fprint("")
        self.fprint(self.conv_rate)
        self.fprint("")

        self.fprint("Taylor Test Finished",special="footer")



        return self.conv_rate

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