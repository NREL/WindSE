from pyadjoint.block import Block
from . import Constant, MPI

class InterpBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x, y, func):
        super(InterpBlock, self).__init__()
        self.x = x
        self.y = y
        self.func = func
        # print(f"{MPI.comm_world.Get_rank()}, {float(x)}, {float(y)}")
        self.add_dependency(x)
        self.add_dependency(y)

    def __str__(self):
        return "InterpBlock"

    def prepare_recompute_component(self, inputs, relevant_outputs):
        x = inputs[0]
        y = inputs[1]
        return [x, y]

    def recompute_component(self, inputs, block_variable, idx, prepared):
        x = prepared[0]
        y = prepared[1]
        # print(f"{MPI.comm_world.Get_rank()}, {float(x)}, {float(y)}, {self.func(x,y)}")

        return Constant(self.func(x,y))

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        x = inputs[0]
        y = inputs[1]
        return [x, y]

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        x = prepared[0]
        y = prepared[1]
        # print(f"{MPI.comm_world.Get_rank()}, {float(x)}, {float(y)}")
        adj_input = adj_inputs[0]

        # Compute derivative with respect to x
        if idx == 0:
            adj = self.func(x,y,dx=1)

        # Compute derivative with respect to y
        elif idx == 1:
            adj = self.func(x,y,dy=1)

        adj_output = adj_input * adj

        return adj_output
