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
    import copy
    from sys import platform
    import time
    from memory_profiler import memory_usage

    ### Import the cumulative parameters ###
    from windse import windse_parameters
    import windse.objective_functions as obj_funcs

    ### Check if we need dolfin_adjoint ###
    if windse_parameters.dolfin_adjoint:
        from dolfin_adjoint import *
        from pyadjoint import tape

    ### This import improves the plotter functionality on Mac ###
    if platform == 'darwin':
        import matplotlib
        matplotlib.use('TKAgg')
    import matplotlib.pyplot as plt
else:
    InequalityConstraint = object
    
import openmdao.api as om


class ObjComp(om.ExplicitComponent):
    """
    OpenMDAO component to wrap the objective computation from dolfin.
    
    Specifically, we use the J and dJ (function and Jacobian) methods
    to compute the function value and derivative values as needed by the
    OpenMDAO optimizers.    
    """
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
    """
    OpenMDAO component to wrap the constraint computation.
    
    A small wrapper used on the fenics methods for computing constraint
    and Jacobian values using the OpenMDAO syntax.
    """
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
    """
    Helper function to gather constraint Jacobians. Adapated from fenics.
    """
    if isinstance(m, list):
        return list(map(gather, m))
    elif hasattr(m, "_ad_to_list"):
        return m._ad_to_list(m)
    else:
        return m  # Assume it is gathered already
    
def om_wrapper(J, initial_DVs, dJ, H, bounds, **kwargs):
    """
    Custom optimization wrapper to use OpenMDAO optimizers with dolfin-adjoint.
    
    Follows the API as defined by dolfin-adjoint.
    
    Parameters
    ----------
    J : object
        Function to compute the model analysis value at a design point.
    initial_DVs : array
        The initial design variables so we can get the array sizing correct
        for the OpenMDAO implementation.
    dJ : object
        Function to compute the Jacobian at a design point.
    H : object
        Function to compute the Hessian at a design point (not used).
    bounds : array
        Array of lower and upper bound values for the design variables.
        
    Returns
    -------
    DVs : array
        The optimal design variable values.
    """
    
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
        folder_output = kwargs["options"]["folder"]
        prob.driver.opt_settings["Summary file"] = os.path.join(folder_output, "SNOPT_summary.out")
        prob.driver.opt_settings["Print file"] = os.path.join(folder_output, "SNOPT_print.out")
        
        if kwargs["options"]["verify_snopt"]:
            prob.driver.opt_settings["Verify level"] = 0
        else:
            prob.driver.opt_settings["Verify level"] = -1
            
        # prob.driver.opt_settings["Major iterations limit"] = 3

    prob.model.add_design_var('DVs', lower=lower_bounds, upper=upper_bounds)
    prob.model.add_objective('obj', ref=kwargs["options"]["obj_ref"], ref0=kwargs["options"]["obj_ref0"])
    
    for idx, constraint_type in enumerate(constraint_types):
        con_name = f'con_{idx}'
        if constraint_type == "eq":
            prob.model.add_constraint(con_name, equals=0.)
        else:
            # Inequality means it's positive from scipy and dolfin
            prob.model.add_constraint(con_name, lower=0.)            


    prob.setup()
    
    prob.set_val('DVs', initial_DVs)
    
    # Optional debugging step to check the total derivatives
    if kwargs["options"]["check_totals"]:
        prob.run_model()
        prob.check_totals()

    else:
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
        self.twist_range = float(self.params["optimization"]["twist_range"])

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

        # self.get_minimum_distance_constraint_func(self.controls, 2*np.mean(self.problem.farm.HH))

        self.fprint("Define Optimizing Functional")
        self.J = self.solver.J
        self.Jhat = ReducedFunctional(self.J, self.controls, eval_cb_post=self.ReducedFunctionalCallback)

        self.Jcurrent = self.J

        self.fprint("Number of Controls: {:d}".format(len(self.controls)),special="header")
        self.OptPrintFunction(self.init_vals,None)
        self.fprint("",special="footer")

        self.SetupConstraints()
        self.fprint("Optimizer Setup",special="footer")

    def SetupConstraints(self):
        ### Build constraints ###
        constraints = []

        ### Pop off the layout specific constraint first ###
        if "min_dist" in self.constraint_types.keys():
            min_dist_dict = self.constraint_types.pop("min_dist")
            if "layout" in self.control_types:
                x_inds = self.indexes[0]
                y_inds = self.indexes[1] 
                RD = np.max(self.farm.get_rotor_diameters())
                constraints.append(MinimumDistanceConstraint(x_inds, y_inds, min_dist_dict["target"]*RD, min_dist_dict["scale"]))

        ### Iterate over remaining objective based constraints
        for key, value in self.constraint_types.items():
            constraints.append(ObjectiveConstraint(self.solver, self.controls, key, value))

        ### Merge constraints into one since pyadjoint has a bug when handling more than one
        if len(constraints) > 0:
            self.merged_constraint = MergedConstraint(constraints,self.fprint)
        else:
            self.merged_constraint = []

        ### Evaluate once
        # self.merged_constraint.function(self.controls)
        # self.merged_constraint.jacobian(self.controls)

    @no_annotations
    def DebugOutput(self,iteration=0, m=[]):
        if self.debug_mode:
            if iteration == 0:
                self.tag_output("n_controls", len(self.controls))
                self.tag_output("obj_value0", float(self.J))

                ### Output initial control values ###
                for i, val in enumerate(self.controls):
                    self.tag_output("val0_"+self.names[i],float(val.values()))

                ### Output gradient ###
                if hasattr(self,"gradients"):
                    for i, d in enumerate(self.gradients):
                        self.tag_output("grad_"+self.names[i],float(d))
                
                ### TODO: Output taylor convergence data
                if hasattr(self,"conv_rate"):
                    pass

            else:

                ### Save the new controls ###
                for k, val in enumerate(m):
                    self.tag_output("val%d_%s" % (self.iteration, self.names[k]),float(val))

                ### Save the objective value ###
                self.tag_output("obj_value%d" % (self.iteration),self.Jcurrent)

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
        self.indexes = [[],[],[],[],[],[],[],[]]
        self.init_vals = []
        j = 0

        if "layout" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[0].append(j)
                j+=1
                self.names.append("x_"+repr(i))
                self.controls.append(Control(self.farm.turbines[i].mx))
                self.control_pointers.append(("x",i))
                self.init_vals.append(Constant(float(self.farm.turbines[i].mx)))

                self.indexes[1].append(j)
                j+=1
                self.names.append("y_"+repr(i))
                self.controls.append(Control(self.farm.turbines[i].my))
                self.control_pointers.append(("y",i))
                self.init_vals.append(Constant(float(self.farm.turbines[i].my)))

        if "yaw" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[2].append(j)
                j+=1
                self.names.append("yaw_"+repr(i))
                self.controls.append(Control(self.farm.turbines[i].myaw))
                self.control_pointers.append(("yaw",i))
                self.init_vals.append(Constant(float(self.farm.turbines[i].myaw)))

        if "axial" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[3].append(j)
                j+=1
                self.names.append("axial_"+repr(i))
                self.controls.append(Control(self.farm.turbines[i].maxial))
                self.control_pointers.append(("axial",i))
                self.init_vals.append(Constant(float(self.farm.turbines[i].maxial)))

        if "lift" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.turbines[0].num_actuator_nodes):
                    self.control_pointers.append(("lift",[i,k]))
                    self.indexes[4].append(j)
                    j+=1
                    self.names.append("lift_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.turbines[i].mlift[k]))
                    self.init_vals.append(Constant(float(self.farm.turbines[i].mlift[k])))

        if "drag" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.turbines[0].num_actuator_nodes):
                    self.control_pointers.append(("drag",[i,k]))
                    self.indexes[5].append(j)
                    j+=1
                    self.names.append("drag_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.turbines[i].mdrag[k]))
                    self.init_vals.append(Constant(float(self.farm.turbines[i].mdrag[k])))

        if "chord" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(0,self.farm.turbines[0].num_actuator_nodes):
                    self.control_pointers.append(("chord",[i,k]))
                    self.indexes[6].append(j)
                    j+=1
                    self.names.append("chord_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.turbines[i].mchord[k]))
                    self.init_vals.append(Constant(float(self.farm.turbines[i].mchord[k])))

        if "twist" in self.control_types:
            for i in self.solver.opt_turb_id:
                # for k in range(self.farm.turbines[0].num_actuator_nodes-1):
                for k in range(0,self.farm.turbines[0].num_actuator_nodes):
                    self.control_pointers.append(("twist",[i,k]))
                    self.indexes[6].append(j)
                    j+=1
                    self.names.append("twist_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.turbines[i].mtwist[k]))
                    self.init_vals.append(Constant(float(self.farm.turbines[i].mtwist[k])))

        if "thrust" in self.control_types:
            for i in self.solver.opt_turb_id:
                # for k in range(self.farm.turbines[0].num_actuator_nodes-1):
                for k in range(0,self.farm.turbines[0].num_actuator_nodes):
                    self.control_pointers.append(("thrust",[i,k]))
                    self.indexes[6].append(j)
                    j+=1
                    self.names.append("thrust_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.turbines[i].mthrust[k]))
                    self.init_vals.append(Constant(float(self.farm.turbines[i].mthrust[k])))

        if "spin" in self.control_types:
            for i in self.solver.opt_turb_id:
                # for k in range(self.farm.turbines[0].num_actuator_nodes-1):
                for k in range(0,self.farm.turbines[0].num_actuator_nodes):
                    self.control_pointers.append(("spin",[i,k]))
                    self.indexes[6].append(j)
                    j+=1
                    self.names.append("spin_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.turbines[i].mspin[k]))
                    self.init_vals.append(Constant(float(self.farm.turbines[i].mspin[k])))

        if "body_force" in self.control_types:
                    self.control_pointers.append(("body_force",[-1,-1]))
                    self.indexes[7].append(j)
                    j+=1
                    self.names.append("body_force")
                    self.controls.append(Control(self.problem.mbody_force))
                    self.init_vals.append(Constant(float(self.problem.mbody_force)))

        self.num_controls = len(self.controls)

    def CreateBounds(self):
        lower_bounds = []
        upper_bounds = []

        if "layout" in self.control_types:
            for i in self.solver.opt_turb_id:
                lower_bounds.append(Constant((self.layout_bounds[0][0])))# + self.farm.radius[i])))
                lower_bounds.append(Constant((self.layout_bounds[1][0])))# + self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[0][1])))# - self.farm.radius[i])))
                upper_bounds.append(Constant((self.layout_bounds[1][1])))# - self.farm.radius[i])))

        if "yaw" in self.control_types:
            for i in self.solver.opt_turb_id:
                lower_bounds.append(Constant(-45*pi/180.0))
                upper_bounds.append(Constant(45*pi/180.0))

        if "axial" in self.control_types:
            for i in rself.solver.opt_turb_id:
                lower_bounds.append(Constant(0))
                upper_bounds.append(Constant(1.))

        if "lift" in self.control_types:
            for i in self.solver.opt_turb_id:
                num_actuator_nodes = self.farm.turbines[i].num_actuator_nodes
                for k in range(num_actuator_nodes):
                    lower_bounds.append(Constant(0))
                    upper_bounds.append(Constant(2.))

        if "drag" in self.control_types:
            for i in self.solver.opt_turb_id:
                num_actuator_nodes = self.farm.turbines[i].num_actuator_nodes
                for k in range(num_actuator_nodes):
                    lower_bounds.append(Constant(0))
                    upper_bounds.append(Constant(2.))

        if "chord" in self.control_types:
            for i in self.solver.opt_turb_id:
                c_avg = 0
                num_actuator_nodes = self.farm.turbines[i].num_actuator_nodes
                max_chord = self.farm.turbines[i].max_chord
                baseline_chord = self.farm.turbines[i].baseline_chord
                for k in range(num_actuator_nodes):
                    modifier = 2.0
                    max_chord = max_chord
                    seg_chord = baseline_chord[k]
                    lower_bounds.append(Constant(seg_chord/modifier))
                    upper_bounds.append(Constant(np.maximum(np.minimum(seg_chord*modifier,max_chord),c_avg)))
                    c_avg = (c_avg*k+seg_chord)/(k+1)

        if "twist" in self.control_types:
            for i in self.solver.opt_turb_id:
                num_actuator_nodes = self.farm.turbines[i].num_actuator_nodes
                twist = self.farm.turbines[i].mtwist
                for k in range(num_actuator_nodes):
                    envelope = np.radians(float(self.twist_range))
                    lower_bounds.append(Constant(float(twist[k])-envelope))
                    upper_bounds.append(Constant(float(twist[k])+envelope))

        if "thrust" in self.control_types:
            for i in self.solver.opt_turb_id:
                num_actuator_nodes = self.farm.turbines[i].num_actuator_nodes
                for k in range(num_actuator_nodes):
                    lower_bounds.append(Constant(0.01))
                    upper_bounds.append(Constant(3.00))

        if "spin" in self.control_types:
            for i in self.solver.opt_turb_id:
                num_actuator_nodes = self.farm.turbines[i].num_actuator_nodes
                for k in range(num_actuator_nodes):
                    lower_bounds.append(Constant(-0.300))
                    upper_bounds.append(Constant(0.300))

        if "body_force" in self.control_types:
            lower_bounds.append(Constant(0.0001))
            upper_bounds.append(Constant(0.01))

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

        self.fprint('========Gradient Data========')
        d_out = '%12s: %12s, %22s' % ('Control', 'Value', 'Derivative')
        self.fprint(d_out)

        d_global = np.zeros(self.params.num_procs, dtype=np.float64)

        for i, d in enumerate(der):
            ctl_val = float(self.controls[i].values())
            # d_out = str(self.names[i] + "    " +repr(ctl_val)+ "    " +repr(float(d)))

            # d_format = np.float64(d)
            # self.params.comm.Gather(d_format, d_global, root=0)
            # d_sum = np.sum(d_global)

            if "yaw" in self.names[i]:
                ctl_val = np.degrees(ctl_val)
                d = d*np.pi/180

            d_out = '%12s: %12.5e, %22.15e' % (self.names[i], ctl_val, d)

            # print('Rank %d, %s' % (self.params.rank, d_out))
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
        # self.fprint("Previous Control Values",special="header")
        # for i, [l, ix] in enumerate(self.control_pointers):
        #     if not isinstance(ix,int):
        #         self.fprint(self.names[i] +": " +repr(float(l[ix[0]][ix[1]])))
        #     else:
        #         self.fprint(self.names[i] +": " +repr(float(l[ix])))
        # self.fprint("",special="footer")

        self.fprint("Next Control Values",special="header")
        for i, val in enumerate(m):
            if "yaw" in self.names[i]:
                val = np.degrees(float(val))
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
        self.problem.farm.save_wind_farm(val=self.iteration)
        # print(m)
        # print(type(m))
        m_new = np.array(m,dtype=float)
        m_old = []
        for control_name, index in self.control_pointers:
            if control_name == "body_force":
                m_old.append(float(self.problem.body_force))
            elif not isinstance(index,int):
                turb_id = index[0]
                seg_id = index[1]
                m_old.append(float(getattr(self.farm.turbines[turb_id],control_name)[seg_id]))
            else:
                m_old.append(float(getattr(self.farm.turbines[index],control_name)))
        # print(m_new)
        # print(type(m_new)) 

        for i in range(len(m_new)):
            if "yaw" in self.names[i]:
                m_new[i] = np.degrees(float(m_new[i]))
                m_old[i] = np.degrees(float(m_old[i]))

        if self.problem.params.rank == 0:
            if self.iteration == 0:

                #### ADD HEADER ####
                self.last_m = np.zeros(self.num_controls)
                for i in range(self.num_controls):
                    self.last_m[i]=float(m_new[i])
                err = 0.0
                f = open(folder_string+"optimization_data.txt",'w')
                header = str("Objective    Change    Prev_Controls:    p_"+"    p_".join(self.names)+"    New_Controls:    n_"+"    n_".join(self.names)+"\n")
                f.write(header)
            else:
                err = np.linalg.norm(m_new-self.last_m)
                self.last_m = copy.copy(m_new)
                f = open(folder_string+"optimization_data.txt",'a')

            output_data = np.concatenate(((self.Jcurrent, err, self.num_controls),m_old))
            output_data = np.concatenate((output_data,(self.num_controls,)))
            output_data = np.concatenate((output_data,m_new))

            np.savetxt(f,[output_data])
            f.close()

    def SaveFunctions(self):
        if self.params["problem"]["type"] == "unsteady":
            u = self.problem.u_k
            p = self.problem.p_k
        else:
            u, p = self.problem.up_k.split()

        # tf = project(self.problem.tf,self.problem.fs.V,solver_type='cg',preconditioner_type="hypre_amg")
        if self.iteration == 0:
            self.velocity_file = self.params.Save(u,"velocity",subfolder="OptSeries/",val=self.iteration)
            self.pressure_file = self.params.Save(p,"pressure",subfolder="OptSeries/",val=self.iteration)
            # self.tf_file = self.params.Save(tf,"tf",subfolder="OptSeries/",val=self.iteration)
        else:
            self.params.Save(u,"velocity",subfolder="OptSeries/",val=self.iteration,file=self.velocity_file)
            self.params.Save(p,"pressure",subfolder="OptSeries/",val=self.iteration,file=self.pressure_file)
            # self.params.Save(tf,"tf",subfolder="OptSeries/",val=self.iteration,file=self.tf_file)
    
    @no_annotations
    def OptPrintFunction(self,m,test=None):
        if test is not None:
            print("Hey, this method actually gives us more info")
        # print(np.array(m,dtype=float))
        # print(np.array(self.control_pointers,dtype=float))
        # print(np.array(self.problem.farm.myaw,dtype=float))
        self.SaveControls(m)
        self.ListControls(m)
        self.SaveFunctions()
        self.solver.EvaluateObjective(opt_iter=self.iteration)

        if "layout" in self.control_types or "yaw" in self.control_types:
            self.problem.farm.plot_farm(filename="wind_farm_step_"+repr(self.iteration),objective_value=self.Jcurrent)

        # if "chord" in self.control_types:
        #     c_lower = np.array(self.bounds[0])[self.indexes[6]] 
        #     c_upper = np.array(self.bounds[1])[self.indexes[6]] 
        #     self.problem.farm.plot_chord(filename="chord_step_"+repr(self.iteration),objective_value=self.Jcurrent,bounds=[c_lower,c_upper])

        self.DebugOutput(self.iteration, m)

        self.iteration += 1

    # def get_minimum_distance_constraint_func(self, m_pos, min_distance=200):
    #     if "layout" in self.control_types and len(self.control_types)==1:
    #         self.dist_constraint = MinimumDistanceConstraint(m_pos, min_distance)
    #     else:
    #         print("minimum distance is supported when only optimizing layout")
    #         self.dist_constraint = None

    def Optimize(self):
        self.fprint("Beginning Optimization",special="header")

        if self.opt_type == "minimize":
            opt_function = minimize
        elif self.opt_type == "maximize":
            opt_function = maximize
        else:
            raise ValueError(f"Unknown optimization type: {self.opt_type}")

        options = {
            "disp" : True,
            "folder" : self.params.folder,
            "gtol": 1.0e-8
            }
        if hasattr(self, 'obj_ref'):
            options["obj_ref"] = self.obj_ref
        if hasattr(self, 'obj_ref0'):
            options["obj_ref0"] = self.obj_ref0
        if hasattr(self, 'verify_snopt'):
            options["verify_snopt"] = self.verify_snopt
        if hasattr(self, 'check_totals'):
            options["check_totals"] = self.check_totals
        
        if self.opt_type == "minimize":
            opt_function = minimize
        elif self.opt_type == "maximize":
            opt_function = maximize
        else:
            raise ValueError(f"Unknown optimization type: {self.opt_type}")

        ### optimize 
        if "SNOPT" in self.opt_routine or "OM_SLSQP" in self.opt_routine:
            m_opt=opt_function(self.Jhat, method="Custom", tol=1.0e-8, options = options, constraints = self.merged_constraint, bounds = self.bounds, callback = self.OptPrintFunction, algorithm=om_wrapper, opt_routine=self.opt_routine)
        else:
            m_opt=opt_function(self.Jhat, method=self.opt_routine, tol=1.0e-8, options = options, constraints = self.merged_constraint, bounds = self.bounds, callback = self.OptPrintFunction)

        self.m_opt = m_opt
        
        if self.num_controls == 1:
            self.m_opt = [self.m_opt]
        self.OptPrintFunction(self.m_opt)
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

        self.fprint("Optimization Finished",special="footer")

        return self.m_opt

    def TaylorTest(self):
        
        self.fprint("Beginning Taylor Test",special="header")

        h = []
        for i,c in enumerate(self.controls):
            h.append(Constant(1.0))
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
    def __init__(self, x_inds, y_inds, target=200, scale=1.0):

        self.x_inds = x_inds
        self.y_inds = y_inds
        self.target = float(target)
        self.scale  = float(scale)
        self.name   = "min_dist"
        # print("In mimimum distance constraint")

    def length(self):
        nconstraints = comb(len(self.x_pos),2.)
        return nconstraints

    def function(self, m):
        ieqcons = []
        np_m = np.array(m,dtype=float)
        x = np_m[self.x_inds]
        y = np_m[self.y_inds]
        n = len(x)

        for i in range(n):
            for j in range(n):
                if j > i:
                    con_value = ((x[i] - x[j])**2 + (y[i] - y[j])**2) - self.target**2
                    ieqcons.append(self.scale*con_value)

        arr = np.array(ieqcons)

        # print("In mimimum distance constraint function eval")
        # print("distances: ", arr)
        self.cur_val = sqrt(np.min(ieqcons)/self.scale+self.target**2)

        numClose = 0
        for i in range(len(arr)):
            if arr[i]<0:
                # print(arr[i]*lengthscale)
                numClose +=1
        if numClose >= 1:
            print("Warning: Number of turbines in violation of spacing constraint: "+repr(numClose))
        return np.array(ieqcons)

    def jacobian(self, m):
        ieqcons = []
        np_m = np.array(m,dtype=float)
        x = np_m[self.x_inds]
        y = np_m[self.y_inds]
        n = len(x)

        for i in range(n):
            for j in range(n):
                if j>i:
                    prime_ieqcons = np.zeros(len(m))

                    prime_ieqcons[2 * i] = self.scale*(2 * (x[i] - x[j]))
                    prime_ieqcons[2 * j] = self.scale*(-2 * (x[i] - x[j]))
                    prime_ieqcons[2 * i + 1] = self.scale*(2 * (y[i] - y[j]))
                    prime_ieqcons[2 * j + 1] = self.scale*(-2 * (y[i] - y[j]))

                    ieqcons.append(prime_ieqcons)
        # print("In mimimum distance constraint Jacobian eval")
        return np.array(ieqcons)

class ObjectiveConstraint(InequalityConstraint):
    def __init__(self, solver, controls, key, value):
        if isinstance(key, int):
            objective_name = value["type"]
        else:
            objective_name = key

        self.name = objective_name

        target = value["target"]
        scale = value["scale"]
        obj_kwargs = value["kwargs"]

        self.target = float(target)
        self.solver = solver
        self.controls = controls
        self.scale = float(scale)
        self.obj_kwargs = obj_kwargs

        ### Precalculate constraint
        angle = self.solver.problem.dom.inflow_angle
        self.J = self.scale*(self.objective(self.solver,angle,**self.obj_kwargs)-self.target)
        self.Jhat = ReducedFunctional(self.J, self.controls)

    def length(self):
        return 1

    def function(self, m):
        # Calculate legacy angle
        angle = self.solver.problem.dom.inflow_angle

        # compute objective and subtract target
        self.cur_val = self.objective(self.solver,angle,**self.obj_kwargs)
        J = self.scale*(self.cur_val-self.target)

        # check if violated
        # print(f"evaluating {self.name}: {self.cur_val} with: {(self.target,self.scale,self.obj_kwargs)}")
        if J < 0:
            print(f"Warning: The {self.name} constraint is violated with a value of {self.cur_val} and a target of {self.target}")

        return J

    def jacobian(self, m):

        dJ = self.Jhat.derivative()

        # print(f"getting gradients of {self.name}")
        # print(np.array(dJ, dtype=float))
        return np.array(dJ, dtype=float)

class MergedConstraint(InequalityConstraint):
    def __init__(self,constraint_list,fprint):
        self.constraint_list = constraint_list
        self.fprint = fprint

    def function(self, m):
        self.fprint("Evaluating Constraints",special="header")
        start = time.time()
        out = []
        for con in self.constraint_list:
            val = con.function(m)
            out = np.append(out, val)
            self.fprint(f"Constraint, {con.name}, return with value: {con.cur_val} and target: {con.target}")
        self.fprint("Completed",special="footer")

        stop = time.time()
        self.fprint("Complete: {:1.2f} s".format(stop-start),special="footer")
        return out

    def jacobian(self, m):
        out = np.empty((0,len(m)))
        for con in self.constraint_list:
            out = np.vstack((out, con.jacobian(m)))
            
        return out
