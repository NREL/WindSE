import numpy as np
from pyadjoint.block import Block
from . import Constant, Point, Cell, Function
from . import MPI, File
import ufl

class MpiEvalBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, u, x, comm=MPI.comm_world, base_eval=None):
        super(MpiEvalBlock, self).__init__()
        self.comm = comm
        self.x0 = x
        self.add_dependency(u)
        self.base_eval = base_eval
        self.coeffs = self.explore(self.x0)
        for c in self.coeffs:
            self.add_dependency(c, no_duplicates=True)

    def __str__(self):
        return "MpiEvalBlock"

    def explore(self, form, coeffs=[]):
        if isinstance(form,Constant):
            if form not in coeffs:
                coeffs.append(form)
        else:
            for part in form.ufl_operands:
                self.explore(part,coeffs)
        return coeffs

    def recompute_x0(self):
        replaced_coeffs = {}
        for block_variable in self.get_dependencies():
            coeff = block_variable.output
            c_rep = block_variable.saved_output
            if coeff in self.coeffs:
                replaced_coeffs[coeff] = c_rep

        x0 = ufl.replace(self.x0, replaced_coeffs)
        return x0

    def prepare_recompute_component(self, inputs, relevant_outputs):
        return None

    def recompute_component(self, inputs, block_variable, idx, prepared):
        u = inputs[0]
        x0 = self.recompute_x0()
        return self.base_eval(u, x0, comm=self.comm)

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        return None

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        x0 = self.recompute_x0()
        p = Point(np.array(x0))
        V = inputs[0].function_space()
        dofs = V.dofmap()
        mesh = V.mesh()
        element = V.element()
        visited = []
        adj_vec = Function(V).vector().get_local()
        adj_vec_size = len(adj_vec)

        adj_input = adj_inputs[idx]

        # print(self.comm.Get_rank(),adj_inputs[idx])
        for cell_idx in range(len(mesh.cells())):
            cell = Cell(mesh, cell_idx)
            if cell.contains(p):
                # print("yes", cell.midpoint().array())
                for ref_dof, dof in enumerate(dofs.cell_dofs(cell_idx)):
                    # print("dof",dof, self.comm.Get_rank(), inputs[0].vector().get_local()[dof])
                    if dof in visited:
                        continue
                    visited.append(dof)
                    basis = element.evaluate_basis(ref_dof,
                                                   p.array(),
                                                   cell.get_coordinate_dofs(),
                                                   cell.orientation())
                    # print(self.comm.Get_rank(),basis)
                    # print(self.comm.Get_rank(),basis.dot(adj_inputs[idx]))
                    if dof < adj_vec_size:
                        adj_vec[dof] = basis.dot(adj_input)



        output = Function(V).vector()
        output[:] = adj_vec

        if any(np.isnan(adj_vec)):
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print("found the nans")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")

            junk=Function(V)
            junk.vector()[:] = adj_vec
            File("mpi_eval.pvd")<<junk

        # print(dir(output))
        # exit()
        # print(f"mpi: {self.comm.Get_rank()}, {output.min()}, {output.max()}")
        # print(self.comm.Get_rank(), min(adj_vec), max(adj_vec))
        # exit()
        return output
