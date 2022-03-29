import numpy as np
from pyadjoint.block import Block
from . import Constant, Point, Cell, Function
from . import MPI

class MpiEvalBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, u, x, comm=MPI.comm_world, base_eval=None):
        super(MpiEvalBlock, self).__init__()
        self.comm = comm
        self.x0 = x
        self.add_dependency(u)
        self.base_eval = base_eval

    def __str__(self):
        return "MpiEvalBlock"

    def prepare_recompute_component(self, inputs, relevant_outputs):
        return None

    def recompute_component(self, inputs, block_variable, idx, prepared):
        u = inputs[0]
        return self.base_eval(u, self.x0, comm=self.comm)

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        return None

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        p = Point(np.array(self.x0))
        V = inputs[0].function_space()
        dofs = V.dofmap()
        mesh = V.mesh()
        element = V.element()
        visited = []
        adj_vec = Function(V).vector().get_local()

        for cell_idx in range(len(mesh.cells())):
            cell = Cell(mesh, cell_idx)
            if cell.contains(p):
                for ref_dof, dof in enumerate(dofs.cell_dofs(cell_idx)):
                    if dof in visited:
                        continue
                    visited.append(dof)
                    basis = element.evaluate_basis(ref_dof,
                                                   p.array(),
                                                   cell.get_coordinate_dofs(),
                                                   cell.orientation())
                    adj_vec[dof] = basis.dot(adj_inputs[idx])

        output = Function(V).vector()
        output[:] = adj_vec

        return output
