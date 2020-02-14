import dolfin
import dolfin_adjoint
import numpy as np
from sys import platform
import os,shutil
import time
import copy
import ufl

from windse.helper_functions import BaseHeight as backend_BaseHeight
from windse.helper_functions import CalculateDiskTurbineForces as backend_CalculateDiskTurbineForces
from windse.helper_functions import CalculateActuatorLineTurbineForces as backend_CalculateActuatorLineTurbineForces
from windse import windse_parameters
from pyadjoint.tape import get_working_tape, annotate_tape, stop_annotating
from pyadjoint.block import Block
from pyadjoint.overloaded_type import create_overloaded_object
from pyadjoint.enlisting import Enlist
import cProfile

### This import improves the plotter functionality on Mac ###
if platform == 'darwin':
    import matplotlib
    matplotlib.use('TKAgg')
import matplotlib.pyplot as plt





def BaseHeight(x,y,ground,**kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)
    with stop_annotating():
        output = backend_BaseHeight(x,y,ground)
    output = create_overloaded_object(output)

    if annotate:        
        block = BaseHeightBlock(x, y, ground, output)
        block.add_output(output.create_block_variable())

        tape = get_working_tape()
        tape.add_block(block)

    return output

class BaseHeightBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x, y, ground, output):
        super(BaseHeightBlock, self).__init__()
        self.x = x
        self.y = y
        self.ground = ground
        self.output = output
        self.add_dependency(x)
        self.add_dependency(y)

    def __str__(self):
        return "BaseHeightBlock"

    def prepare_recompute_component(self, inputs, relevant_outputs):
        x = inputs[0]
        y = inputs[1]
        return [x, y]

    def recompute_component(self, inputs, block_variable, idx, prepared):
        x = prepared[0]
        y = prepared[1]
        return backend_BaseHeight(x,y,self.ground)

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        x = inputs[0]
        y = inputs[1]
        return [x, y]

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        x = prepared[0]
        y = prepared[1]
        # h = dolfin.Constant(10)
        adj_input = adj_inputs[0]

        # print((float(x),float(y)))
        # print(self.ground(x,y))
        # print((self.ground(x,y,dx=1),self.ground(x,y,dy=1)))

        if idx == 0:
            # print("taking x derivative")
            adj = self.ground(x,y,dx=1)
            # adj = (self.ground(x+h,y)-self.ground(x-h,y))/(2*h)
        elif idx == 1:
            # print("taking y derivative")
            adj = self.ground(x,y,dy=1)
            # adj = (self.ground(x,y+h)-self.ground(x,y-h))/(2*h)


        return adj_input * adj














def CalculateDiskTurbineForces(x,farm,fs,**kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)

    ### Get the actual output ###
    with stop_annotating():
        func_output, sparse_ids, actuator_array = backend_CalculateDiskTurbineForces(x,farm,fs,**kwargs)

    ### Create the block object ###
    if annotate:
        block = ActuatorDiskForceBlock(x,farm,fs,sparse_ids,**kwargs)

    ### Prepare outputs for recording on tape ###
    r = []
    func_output = Enlist(func_output)
    for out in func_output:
        if not isinstance(out,list):
            r.append(create_overloaded_object(out))

    ### Add everything to tape ###
    if annotate:
        for out in r:
            if not isinstance(out,list):
                block.add_output(out.create_block_variable())
        tape = get_working_tape()
        tape.add_block(block)

    return func_output.delist(r), sparse_ids, actuator_array

class ActuatorDiskForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x,farm,fs,sparse_ids,**kwargs):
        super(ActuatorDiskForceBlock, self).__init__()
        self.x = x
        self.farm = farm
        self.fs = fs
        self.sparse_ids = sparse_ids
        self.inflow_angle = copy.copy(farm.inflow_angle)
        self.save_actuator = kwargs.get('save_actuator',False)
        self.sparse_RDs = kwargs.get('sparse_RDs',2.0)
        self.dim = farm.dom.dim

        ### Calculate how far the turbines need to move to recalculate sparse ids ###
        self.move_tol = 0.1*farm.dom.mesh.rmin()

        ### Add dependencies on the controls ###
        for i in range(self.farm.numturbs):
            self.farm.mx[i].block_variable.tag = ("x",i)
            self.add_dependency(self.farm.mx[i])

            self.farm.my[i].block_variable.tag = ("y",i)
            self.add_dependency(self.farm.my[i])

            self.farm.myaw[i].block_variable.tag = ("yaw",i)
            self.add_dependency(self.farm.myaw[i])

            self.farm.ma[i].block_variable.tag = ("a",i)
            self.add_dependency(self.farm.ma[i])

        ### Tabulate which controls matter ###
        self.control_types = []
        if "layout" in farm.control_types:
            self.control_types.append("x")
            self.control_types.append("y")
        if "yaw" in farm.control_types:
            self.control_types.append("yaw")
        if "axial" in farm.control_types:
            self.control_types.append("a")

        ### Record the current values for x and y ###
        self.old_x = np.array(self.farm.mx,dtype=float)
        self.old_y = np.array(self.farm.my,dtype=float)

        ### map ids to vector ids ###
        self.VectorSparseIDs()

    def __str__(self):
        return "ActuatorLineForceBlock"

    def CheckTurbineLocations(self,x,y):
        new_x = np.array(x,dtype=float)
        new_y = np.array(y,dtype=float)
        moved = max(max(abs(new_x-self.old_x)),max(abs(new_y-self.old_y)))
        # print(moved,self.move_tol, moved > self.move_tol)
        if moved > self.move_tol:
            self.sparse_ids = None
        self.old_x = new_x
        self.old_y = new_y

    def VectorSparseIDs(self):
        self.all_ids = np.zeros(len(self.sparse_ids)*self.dim,dtype=int)
        self.all_ids[0::self.dim] = self.dim*(self.sparse_ids)+0
        self.all_ids[1::self.dim] = self.dim*(self.sparse_ids)+1
        if self.dim == 3:
            self.all_ids[2::self.dim] = self.dim*(self.sparse_ids)+2
        self.all_ids = list(self.all_ids)

    def prepare_recompute_component(self, inputs, relevant_outputs):
        ### update the new controls inside the windfarm ###
        x = inputs[0::4]
        y = inputs[1::4]
        yaw = inputs[2::4]
        a = inputs[3::4]
        self.farm.inflow_angle = self.inflow_angle
        self.farm.UpdateControls(x=x,y=y,yaw=yaw,a=a)

        ### determine if we need to recalculate sparse_ids ###
        self.CheckTurbineLocations(x,y)

        prepared, self.sparse_ids, actuator_array = backend_CalculateDiskTurbineForces(self.x,self.farm,self.fs,sparse_ids=self.sparse_ids,sparse_RDs=self.sparse_RDs)

        return prepared

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared[idx]

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        ### update the new controls inside the windfarm ###
        x = inputs[0::4]
        y = inputs[1::4]
        yaw = inputs[2::4]
        a = inputs[3::4]
        self.farm.inflow_angle = self.inflow_angle
        self.farm.UpdateControls(x=x,y=y,yaw=yaw,a=a)

        ### determine if we need to recalculate sparse_ids ###
        self.CheckTurbineLocations(x,y)

        prepared = {}
        for name in self.control_types:
            # pr = cProfile.Profile()
            # pr.enable()
            # print("calculating "+name+" derivatives")
            prepared[name], self.sparse_ids, actuator_array = backend_CalculateDiskTurbineForces(self.x,self.farm,self.fs,dfd=name,sparse_ids=self.sparse_ids,sparse_RDs=self.sparse_RDs)
            # pr.disable()
            # pr.print_stats()

            # exit()

        ### map ids to vector ids ###
        self.VectorSparseIDs()

        return prepared

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, turb_idx = block_variable.tag
        print("Calculating Derivative: " +name+"_"+repr(turb_idx))
        ### Apply derivative to previous in tape ###
        adj_output = 0
        for i in range(3):
            adj_output += np.inner(adj_inputs[i].get_local(self.all_ids),prepared[name][i][:,turb_idx])

        return np.array(adj_output)
















# ================================================================
# ================================================================

def CalculateActuatorLineTurbineForces(problem, simTime, **kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)

    # Example of how to call the update actuator force function
    # UpdateActuatorLineForce(problem, simTime, dfd, tf):

    ### Get the actual output ###
    with stop_annotating():
        func_output = backend_CalculateActuatorLineTurbineForces(problem, simTime, **kwargs)

    ### Create the block object ###
    if annotate:
        block = ActuatorLineForceBlock(problem, simTime,**kwargs)

    ### Prepare outputs for recording on tape ###
    r = []
    func_output = Enlist(func_output)
    for out in func_output:
        if not isinstance(out,list):
            r.append(create_overloaded_object(out))

    ### Add everything to tape ###
    if annotate:
        for out in r:
            if not isinstance(out,list):
                block.add_output(out.create_block_variable())
        tape = get_working_tape()
        tape.add_block(block)

    return func_output.delist(r)

# ================================================================
# ================================================================

class ActuatorLineForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, problem, simTime, **kwargs):
        super(ActuatorLineForceBlock, self).__init__()
        self.problem = problem
        self.simTime = simTime

        # Add dependencies on the controls
        for i in range(self.problem.num_blade_segments):
            self.problem.mcl[i].block_variable.tag = ("c_lift", i)
            self.add_dependency(self.problem.mcl[i])

            self.problem.mcd[i].block_variable.tag = ("c_drag", i)
            self.add_dependency(self.problem.mcd[i])

        # Tabulate which controls matter
        self.control_types = []
        if "lift" in problem.farm.control_types:
            self.control_types.append("c_lift")
        if "drag" in problem.farm.control_types:
            self.control_types.append("c_drag")

        # Record the current values for controls
        self.old_c_lift = np.array(self.problem.mcl, dtype=float)
        self.old_c_drag = np.array(self.problem.mcd, dtype=float)


    def __str__(self):
        return "ActuatorLineForceBlock"


    def prepare_recompute_component(self, inputs, relevant_outputs):
        # update the new controls inside the windfarm
        c_lift = inputs[0::2]
        c_drag = inputs[1::2]

        self.problem.UpdateActuatorLineControls(c_lift = c_lift, c_drag = c_drag)

        prepared = backend_CalculateActuatorLineTurbineForces(self.problem, self.simTime)

        return prepared

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        ### update the new controls inside the windfarm ###
        c_lift = inputs[0::2]
        c_drag = inputs[1::2]

        self.problem.UpdateActuatorLineControls(c_lift = c_lift, c_drag = c_drag)

        prepared = {}
        for name in self.control_types:
            # pr = cProfile.Profile()
            # pr.enable()
            # print("calculating "+name+" derivatives")
            prepared[name] = backend_CalculateActuatorLineTurbineForces(self.problem, self.simTime, dfd = name)
            # pr.disable()
            # pr.print_stats()

            # exit()

        return prepared

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, segment_index = block_variable.tag
        # print("calculating: " + name + "_" + repr(segment_index))

        ### Get the vectors ###
        comp_vec = prepared[name][:, segment_index]
        adj_vec = adj_inputs[0].get_local()

        ### Apply derivative to previous in tape ###
        adj_output = np.inner(comp_vec, adj_vec)
        return np.array(adj_output)

# ================================================================
# ================================================================
























#################################################################################################
############################## Overridden Dolfin_Adjoint Functions ##############################
#################################################################################################


### We need to override ALE move so that dolfin-adjoint can track the change ###
backend_move = dolfin.ALE.move
def move(mesh,bmesh, **kwargs):
    """ Refine is overloaded to ensure that the returned mesh is overloaded.
    """
    with stop_annotating():
        backend_move(mesh,bmesh, **kwargs)
    overloaded = create_overloaded_object(mesh)
    return overloaded
dolfin.ALE.move = move

def linalg_solve(*args, **kwargs):
    """This function overrides dolfin_adjoints.compat.linalg_solve.

    The original function doesn't allow for solver options because it uses
     the::

        dolfin.solve(A,x,b) 

    form which doesn't accept keyword arguments. However, It does except 
    additional arguments that defined some solver options, which we pass
    in manually

    Todo:

        Eventually, we want to replace this with a full PetscKrylovSolver()
        to get access to all the ksp options.

    """
    print("performing a solve")
    return dolfin_adjoint.backend.solve(*args)#,"mumps") 
dolfin_adjoint.types.compat.linalg_solve = linalg_solve

def assemble_adjoint_value(form, **kwargs):
    """Wrapper that assembles a matrix with boundary conditions"""
    bcs = kwargs.pop("bcs", ())
    # print(form)
    if windse_parameters["wind_farm"].get("turbine_method","dolfin") == "numpy":
        rep = 'tsfc'
    else:
        rep = 'uflacs'
    result = dolfin_adjoint.backend.assemble(form,form_compiler_parameters={'representation': rep})
    for bc in bcs:
        bc.apply(result)
    return result
dolfin_adjoint.types.compat.assemble_adjoint_value = assemble_adjoint_value

shutil.rmtree(windse_parameters.folder+"debug/", ignore_errors=True)
def recompute_component(self, inputs, block_variable, idx, prepared):
    file_exists = False
    # print("resolving")

    # print()
    # print("hey Look at me I'm recomputing!")
    # print()
    """This function overrides 
    dolfin_adjoint.solving.SolveBlock.recompute_component

    The original function doesn't account for kwargs in a project so it 
    doesn't pass them to this function. For now, we just supply the
    solver_parameters manually in this case. 

    Todo:

        Eventually, we want to replace this with a full PetscKrylovSolver()
        to get access to all the ksp options.

    """
    lhs = prepared[0]
    rhs = prepared[1]
    func = prepared[2]
    bcs = prepared[3]


    if not self.forward_kwargs:
        dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, solver_parameters={'linear_solver': 'mumps'})
    else:
        dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, **self.forward_kwargs)

        if "debug" in windse_parameters.output:
            if hasattr(self, 'solve_iteration'):
                self.solve_iteration += 1
            else:
                # print()
                # print("but we haven't been here before")
                self.recompute_set = 0
                while os.path.isfile(windse_parameters.folder+"debug/dolfin_adjoint_func_"+repr(self.recompute_set)+".pvd"):
                    self.recompute_set +=1
                # print("and that's the "+repr(self.recompute_set)+" time this has happened")
                # print()
                self.savefile = dolfin.File(windse_parameters.folder+"debug/dolfin_adjoint_func_"+repr(self.recompute_set)+".pvd")
                self.solve_iteration = 0
            
            
            u, p = func.split(True)
            u.rename("velocity","velocity")
            self.savefile << (u,self.solve_iteration)
    # print("assemble(func*dx): " + repr(float(dolfin.assemble(dolfin.inner(func,func)*dolfin.dx))))
    return func
dolfin_adjoint.solving.SolveBlock.recompute_component = recompute_component

def _init_dependencies(self, *args, **kwargs):
    self.preconditioner_method = "default"
    self.solver_method = "mumps"

    if self.varform:
        eq = args[0]
        self.lhs = eq.lhs
        self.rhs = eq.rhs
        self.func = args[1]

        if len(args) > 2:
            self.bcs = args[2]
        elif "bcs" in kwargs:
            self.bcs = self.kwargs.pop("bcs")
            self.forward_kwargs.pop("bcs")
        else:
            self.bcs = []

        if self.bcs is None:
            self.bcs = []

        self.assemble_system = False
    else:
        # Linear algebra problem.
        # TODO: Consider checking if attributes exist.
        A = args[0]
        u = args[1]
        b = args[2]
        if len(args) >= 4:
            self.solver_method = args[3]
        if len(args) >= 5:
            self.preconditioner_method = args[4]

        sp = {"solver_parameters": {'linear_solver': self.solver_method,
                                    'preconditioner': self.preconditioner_method}}
        self.forward_kwargs.update(sp)

        self.lhs = A.form
        self.rhs = b.form
        self.bcs = A.bcs if hasattr(A, "bcs") else []
        self.func = u.function
        self.assemble_system = A.assemble_system if hasattr(A, "assemble_system") else False

    if not isinstance(self.bcs, list):
        self.bcs = [self.bcs]

    if isinstance(self.lhs, ufl.Form) and isinstance(self.rhs, ufl.Form):
        self.linear = True
        # Add dependence on coefficients on the right hand side.
        for c in self.rhs.coefficients():
            self.add_dependency(c, no_duplicates=True)
    else:
        self.linear = False

    for bc in self.bcs:
        self.add_dependency(bc, no_duplicates=True)

    for c in self.lhs.coefficients():
        self.add_dependency(c, no_duplicates=True)
dolfin_adjoint.solving.SolveBlock._init_dependencies = _init_dependencies

def _assemble_and_solve_adj_eq(self, dFdu_form, dJdu):
    dJdu_copy = dJdu.copy()
    kwargs = self.assemble_kwargs.copy()
    # Homogenize and apply boundary conditions on adj_dFdu and dJdu.
    bcs = self._homogenize_bcs()
    kwargs["bcs"] = bcs
    dFdu = dolfin_adjoint.compat.assemble_adjoint_value(dFdu_form, **kwargs)

    for bc in bcs:
        bc.apply(dJdu)

    adj_sol = dolfin_adjoint.Function(self.function_space)
    dolfin_adjoint.compat.linalg_solve(dFdu, adj_sol.vector(), dJdu, self.solver_method, self.preconditioner_method,**self.kwargs)

    adj_sol_bdy = dolfin_adjoint.compat.function_from_vector(self.function_space, dJdu_copy - dolfin_adjoint.compat.assemble_adjoint_value(
        dolfin.action(dFdu_form, adj_sol)))

    return adj_sol, adj_sol_bdy
dolfin_adjoint.solving.SolveBlock._assemble_and_solve_adj_eq = _assemble_and_solve_adj_eq








