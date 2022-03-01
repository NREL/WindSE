from pyadjoint.block import Block

class MarkerBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, u, simTime, out_file):
        super(MarkerBlock, self).__init__()
        self.simTime = simTime
        self.out_file = out_file
        self.saved = False
        self.V = u.function_space() 
        self.add_dependency(u)

    def __str__(self):
        return "MarkerBlock"

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return inputs[0]

    # @no_annotations
    # def evaluate_adj(self, markings=False):
    #     print("Saving Adjoint Function for timestep:" + repr(self.simTime))

    #     output = self.get_outputs()
    #     adj_input = output[0].adj_value

    #     temp = dolfin_adjoint.Function(self.V, name="u_adjoint",annotate=False)
    #     temp.vector().set_local(adj_input.get_local())
    #     self.out_file.write(temp,self.simTime)



    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        if self.saved == False:
            print("Saving Adjoint Function for timestep: " + repr(self.simTime))
            temp = dolfin_adjoint.Function(self.V, name="u_adjoint",annotate=False)
            temp.vector().set_local(adj_inputs[0].get_local())
            self.out_file.write(temp,self.simTime)
            self.saved = True

        return adj_inputs[0]

