import dolfin
import dolfin_adjoint
import numpy as np

from windse.helper_functions import BaseHeight as backend_BaseHeight
from pyadjoint.tape import get_working_tape, annotate_tape, stop_annotating
from pyadjoint.block import Block
from pyadjoint.overloaded_type import create_overloaded_object

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
def recompute_component(self, inputs, block_variable, idx, prepared):
    """This function overrides
    dolfin_adjoint.solving.SolveBlock.recompute_component

    The original function doesn't account for kwargs in a project so it
    doesn't pass them to this function. For now, we just supply the
    solver_parameters manually in this case.

    Todo:

        Eventually, we want to replace this with a full PetscKrylovSolver()
        to get access to all the ksp options.

    """
    eq = prepared[0]
    func = prepared[1]
    bcs = prepared[2]


    if not self.forward_kwargs:
        dolfin_adjoint.backend.solve(eq, func, bcs, solver_parameters={'linear_solver': 'mumps'})
    else:
        dolfin_adjoint.backend.solve(eq, func, bcs, **self.forward_kwargs)
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

        tape = get_working_tape()
        tape.add_block(block)

        block.add_output(output.block_variable)

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
        self.add_dependency(x.block_variable)
        self.add_dependency(y.block_variable)

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
        h = dolfin.Constant(1e-5)
        adj_input = adj_inputs[0]

        if idx == 0:
            adj = (self.ground(x+h,y)-self.ground(x-h,y))/(2*h)
        elif idx == 1:
            adj = (self.ground(x,y+h)-self.ground(x,y-h))/(2*h)

        return adj_input * adj

class ReducedFunctional(dolfin_adjoint.ReducedFunctional):
    """This is a special version of the ReducedFunctional class allowing
    windse to inject special function calls.

    Args:
        values ([OverloadedType]): If you have multiple controls this
            should be a list of new values for each control in the order
            you listed the controls to the constructor. If you have a
            single control it can either be a list or a single object.
            Each new value should have the same type as the
            corresponding control.

    Returns:
        :obj:`OverloadedType`: The computed value. Typically of instance
            of :class:`AdjFloat`.

    """
    def __init__(self, *args, **kwargs):
        super(ReducedFunctional, self).__init__(*args, **kwargs)
        self.iter_complete = False

    def __call__(self, values):
        if self.iter_complete:
            self.display_output()
            self.iter_complete = False


        values = dolfin_adjoint.pyadjoint.enlisting.Enlist(values)
        if len(values) != len(self.controls):
            raise ValueError("values should be a list of same length as controls.")

        # Call callback.
        self.eval_cb_pre(self.controls.delist(values))

        for i, value in enumerate(values):
            self.controls[i].update(value)

        blocks = self.tape.get_blocks()
        with self.marked_controls():
            with dolfin_adjoint.pyadjoint.tape.stop_annotating():
                for i in range(len(blocks)):
                    blocks[i].recompute()

        func_value = self.functional.block_variable.checkpoint

        # Call callback
        self.eval_cb_post(func_value, self.controls.delist(values))

        return func_value

    def display_output(self):
        print("$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$")
        print("Custom Code Here!")
        print("$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$")
        print("$$$$$$$$$$$$$$$$$")

















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
