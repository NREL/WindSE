import dolfin
import dolfin_adjoint
import numpy as np
from sys import platform
import os,shutil
import time
import copy
import ufl

from windse import windse_parameters
from pyadjoint.tape import get_working_tape, annotate_tape, stop_annotating, no_annotations
from pyadjoint.block import Block
from pyadjoint.overloaded_type import create_overloaded_object
from pyadjoint.enlisting import Enlist
import cProfile


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

if "dolfin_adjoint_helper" not in dolfin.ALE.move.__module__:
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
    # return dolfin_adjoint.backend.solve(*args,"mumps") 

if "dolfin_adjoint_helper" not in dolfin_adjoint.types.compat.linalg_solve.__module__:
    dolfin_adjoint.types.compat.linalg_solve = linalg_solve

def assemble_adjoint_value(form, **kwargs):
    """Wrapper that assembles a matrix with boundary conditions"""
    bcs = kwargs.pop("bcs", ())
    # print(form)
    if windse_parameters["turbines"]["type"] == "numpy_disk":
        rep = 'tsfc'
    else:
        rep = 'uflacs'
    result = dolfin_adjoint.backend.assemble(form,form_compiler_parameters={'representation': rep})
    for bc in bcs:
        bc.apply(result)
    return result

if "dolfin_adjoint_helper" not in dolfin_adjoint.types.compat.assemble_adjoint_value.__module__:
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
        # dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, solver_parameters={'linear_solver': 'mumps'})
        dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, solver_parameters={'linear_solver': 'gmres', 'preconditioner': 'hypre_amg'})
    else:
        if "krylov_solver_parameters" in self.forward_kwargs.keys():
            solver_parameters={'linear_solver': self.forward_kwargs["krylov_method"], 'preconditioner': self.forward_kwargs["krylov_preconditioner"]}
            dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, solver_parameters=solver_parameters)
        else:
            dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, **self.forward_kwargs)
        # print()
        # print("func = "+repr(np.mean(func.vector().get_local())))

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
if "dolfin_adjoint_helper" not in dolfin_adjoint.solving.SolveBlock.recompute_component.__module__:
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
if "dolfin_adjoint_helper" not in dolfin_adjoint.solving.SolveBlock._init_dependencies.__module__:
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
if "dolfin_adjoint_helper" not in dolfin_adjoint.solving.SolveBlock._assemble_and_solve_adj_eq.__module__:
    dolfin_adjoint.solving.SolveBlock._assemble_and_solve_adj_eq = _assemble_and_solve_adj_eq

# This is a fix to account for a strange occurrence involving as_vector(np_a_float)
def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
    adj_output = np.zeros(inputs[0].shape)
    adj_input = adj_inputs[0]
    if isinstance(adj_input,dolfin.Vector):
        adj_output[self.item] = adj_input.get_local()
    else:
        adj_output[self.item] = adj_input
    return adj_output

if "dolfin_adjoint_helper" not in dolfin_adjoint.numpy_adjoint.array.NumpyArraySliceBlock.evaluate_adj_component.__module__:
    dolfin_adjoint.numpy_adjoint.array.NumpyArraySliceBlock.evaluate_adj_component = evaluate_adj_component

def dot(self, o, da, db):
    a, b = o.ufl_operands
    pa = dolfin.dot(da, b)
    pb = dolfin.dot(a, db)
    s = ufl.classes.Sum(pa, pb)
    return s

if not hasattr(ufl.algorithms.apply_derivatives.GenericDerivativeRuleset, "dot"):
    ufl.algorithms.apply_derivatives.GenericDerivativeRuleset.dot = dot

def cross(self, o, da, db):
    a, b = o.ufl_operands
    pa = dolfin.cross(da, b)
    pb = dolfin.cross(a, db)
    s = ufl.classes.Sum(pa, pb)
    return s

if not hasattr(ufl.algorithms.apply_derivatives.GenericDerivativeRuleset, "cross"):
    ufl.algorithms.apply_derivatives.GenericDerivativeRuleset.cross = cross

def transposed(self, o, da):
    return da.T

if not hasattr(ufl.algorithms.apply_derivatives.GenericDerivativeRuleset, "transposed"):
    ufl.algorithms.apply_derivatives.GenericDerivativeRuleset.transposed = transposed
