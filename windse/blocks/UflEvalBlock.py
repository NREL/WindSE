import numpy as np
from pyadjoint.block import Block
from . import Constant, diff, create_overloaded_object
import ufl

class UflEvalBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, form, print_statement=None, base_eval=None):
        super(UflEvalBlock, self).__init__()
        self.form = form
        self.base_eval = base_eval
        self.coeffs = self.explore(self.form)
        self.print_statement = print_statement
        for c in self.coeffs:
            self.add_dependency(c, no_duplicates=True)

    def __str__(self):
        return "UflEvalBlock"

    def explore(self, form, coeffs=[]):
        if isinstance(form,Constant):
            if form not in coeffs:
                coeffs.append(form)
        else:
            for part in form.ufl_operands:
                self.explore(part,coeffs)
        return coeffs

    def prepare_recompute_component(self, inputs, relevant_outputs):
        return self.prepare_evaluate_adj(inputs, None, None)

    def recompute_component(self, inputs, block_variable, idx, prepared):
        if self.print_statement is not None:
            print(f"Forward: {self.print_statement}")
        form = prepared
        output = self.base_eval(form)
        output = create_overloaded_object(output)
        return output

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        replaced_coeffs = {}
        for block_variable in self.get_dependencies():
            coeff = block_variable.output
            c_rep = block_variable.saved_output
            if coeff in self.coeffs:
                replaced_coeffs[coeff] = c_rep

        form = ufl.replace(self.form, replaced_coeffs)
        return form

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        print("adj_test")
        if self.print_statement is not None:
            print(f"Adjoint: {self.print_statement}")
        form = prepared
        adj_input = adj_inputs[0]
        der = diff(form,inputs[idx])
        print(der)
        output = float(ufl.algorithms.apply_derivatives.apply_derivatives(der))
        print(output)
        return adj_input * output
