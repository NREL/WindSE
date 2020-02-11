import dolfin
import dolfin_adjoint
import numpy as np
from sys import platform
import os,shutil
import time
import copy

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
    return dolfin_adjoint.backend.solve(*args,"mumps") 
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

        # print(eq.lhs)

        # exit()
        # print()
        # print("Oh man and I'm doing a fancy solve")
        # print()
        # test = dolfin.File("tf_x.pvd")
        # test << eq.lhs.coefficients()[3]
        # test = dolfin.File("tf_y.pvd")
        # test << eq.lhs.coefficients()[4]
        # test = dolfin.File("u0.pvd")
        # test << eq.lhs.coefficients()[5]

        # bc_dofs = bcs[0].get_boundary_values().keys()
        # bc_x = []
        # bc_y = []
        # for dof in bc_dofs:
        #     bc_x.append(func.function_space().tabulate_dof_coordinates()[dof][0])
        #     bc_y.append(func.function_space().tabulate_dof_coordinates()[dof][1])

        # print("Inflow Velocity: " + repr(bcs[0].value()([0.0,0.0])))
        # print("Initial Condition: " + repr(func([0.0,0.0])))
        # plt.clf()
        # plt.scatter(bc_x,bc_y)
        # plt.savefig("Marker.pdf")
        # plt.show()

        # exit()

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

    return func

dolfin_adjoint.solving.SolveBlock.recompute_component = recompute_component

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
        block = TurbineForceBlock(x,farm,fs,sparse_ids,**kwargs)

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

class TurbineForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x,farm,fs,sparse_ids,**kwargs):
        super(TurbineForceBlock, self).__init__()
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
        return "TurbineForceBlock"

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

        ### Apply derivative to previous in tape ###
        adj_output = 0
        for i in range(3):
            adj_output += np.inner(adj_inputs[i].get_local(self.all_ids),prepared[name][i][:,turb_idx])

        return np.array(adj_output)



















# class ReducedFunctional(dolfin_adjoint.ReducedFunctional):
#     """This is a special version of the ReducedFunctional class allowing
#     windse to inject special function calls.

#     Args:
#         values ([OverloadedType]): If you have multiple controls this 
#             should be a list of new values for each control in the order
#             you listed the controls to the constructor. If you have a 
#             single control it can either be a list or a single object. 
#             Each new value should have the same type as the 
#             corresponding control.

#     Returns:
#         :obj:`OverloadedType`: The computed value. Typically of instance
#             of :class:`AdjFloat`.

#     """
#     def __init__(self, *args, **kwargs):
#         super(ReducedFunctional, self).__init__(*args, **kwargs)
#         self.iter_complete = False

#     def __call__(self, values):
#         if self.iter_complete:
#             self.display_output()
#             self.iter_complete = False


#         values = dolfin_adjoint.pyadjoint.enlisting.Enlist(values)
#         if len(values) != len(self.controls):
#             raise ValueError("values should be a list of same length as controls.")

#         # Call callback.
#         self.eval_cb_pre(self.controls.delist(values))

#         for i, value in enumerate(values):
#             self.controls[i].update(value)

#         blocks = self.tape.get_blocks()
#         with self.marked_controls():
#             with dolfin_adjoint.pyadjoint.tape.stop_annotating():
#                 for i in range(len(blocks)):
#                     blocks[i].recompute()

#         func_value = self.functional.block_variable.checkpoint

#         # Call callback
#         self.eval_cb_post(func_value, self.controls.delist(values))

#         return func_value

#     def display_output(self):
#         print("$$$$$$$$$$$$$$$$$")
#         print("$$$$$$$$$$$$$$$$$")
#         print("$$$$$$$$$$$$$$$$$")
#         print(dir(self))
#         print("$$$$$$$$$$$$$$$$$")
#         print("$$$$$$$$$$$$$$$$$")
#         print("$$$$$$$$$$$$$$$$$")

















#### This is legacy code that hacked into dolfin_adjoint and forced the base heights to update ####
# def update_relative_heights(tape):
#     """This function find the the turbine (x,y,z) control values and updates
#     the z values according to the updated (x,y) values that are being
#     optimized. 

#     """

#     ### This gets the list of Constants associated with the turbine force ###
#     blocks = tape.get_blocks()
#     depends = blocks[0].get_dependencies()

#     ### This loops over the Constants and identifies them ###
#     x_ind = []
#     y_ind = []
#     z_ind = []
#     for i in range(len(depends)):
#         cur = depends[i]

#         # print(cur.output)
#         # print(cur.saved_output)

#         if "x" in str(cur.output):
#             x_ind.append(i)
#             # x_ind.append(str(cur.saved_output))
#         if "y" in str(cur.output):
#             y_ind.append(i)
#             # y_ind.append(str(cur.saved_output))
#         if "z" in str(cur.output):
#             z_ind.append(i)
#             # z_ind.append(str(cur.saved_output))

#     # print(x_ind)
#     # print(y_ind)
#     # print(z_ind)
#     # print()
#     # print()

#     ### Finally we extract the x and y values and update the z values ###
#     for i in range(len(z_ind)):
#         x_val = depends[x_ind[i]].saved_output.values()[0]
#         y_val = depends[y_ind[i]].saved_output.values()[0]
#         z_val = float(windse_parameters.ground_fx(x_val,y_val)) + windse_parameters.full_hh[i]
#         depends[z_ind[i]].saved_output.assign(z_val)
#         print(x_val,y_val,z_val)

#     # exit()


# def reduced_functional_eval(self, values):
#     """This function overrides 
#     pyadjoint.reduced_functional.ReducedFunctional.__call__() allowing
#     windse to update the turbine absolute heights after the locations
#     are updated.

#     Args:
#         values ([OverloadedType]): If you have multiple controls this 
#             should be a list of new values for each control in the order
#             you listed the controls to the constructor. If you have a 
#             single control it can either be a list or a single object. 
#             Each new value should have the same type as the 
#             corresponding control.

#     Returns:
#         :obj:`OverloadedType`: The computed value. Typically of instance
#             of :class:`AdjFloat`.

#     """
#     values = dolfin_adjoint.pyadjoint.enlisting.Enlist(values)
#     if len(values) != len(self.controls):
#         raise ValueError("values should be a list of same length as controls.")

#     # Call callback.
#     self.eval_cb_pre(self.controls.delist(values))

#     for i, value in enumerate(values):
#         self.controls[i].update(value)

#     ### This is the new code injected into pyadjoint ###
#     update_relative_heights(self.tape)
#     ####################################################

#     blocks = self.tape.get_blocks()
#     with self.marked_controls():
#         with dolfin_adjoint.pyadjoint.tape.stop_annotating():
#             for i in range(len(blocks)):
#                 blocks[i].recompute()

#     func_value = self.functional.block_variable.checkpoint

#     # Call callback
#     self.eval_cb_post(func_value, self.controls.delist(values))

#     return func_value

# dolfin_adjoint.ReducedFunctional.__call__ = reduced_functional_eval