from pyadjoint.block import Block
import numpy as np

class ActuatorDiskForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x, inflow_angle, fs, sparse_ids, turb=None, construct_sparse_ids=None):
        super(ActuatorDiskForceBlock, self).__init__()
        self.x = x
        self.fs = fs
        self.inflow_angle = inflow_angle
        self.sparse_ids = sparse_ids
        self.turb = turb
        self.construct_sparse_ids = construct_sparse_ids
        self.sparse_RDs = 2.0
        self.dim = turb.dom.dim

        ### Calculate how far the turbines need to move to recalculate sparse ids ###
        self.move_tol = 0.1*turb.dom.mesh.rmin()

        ### Add dependencies on the controls ###
        self.turb.mx.block_variable.tag = "x"
        self.add_dependency(self.turb.mx)

        self.turb.my.block_variable.tag = "y"
        self.add_dependency(self.turb.my)

        self.turb.myaw.block_variable.tag = "yaw"
        self.add_dependency(self.turb.myaw)

        self.turb.maxial.block_variable.tag = "axial"
        self.add_dependency(self.turb.maxial)

        ### Record the current values for x and y ###
        self.old_x = self.turb.x
        self.old_y = self.turb.x

        ### map ids to vector ids ###
        self.vector_sparse_ids()

    def __str__(self):
        return "ActuatorLineForceBlock"

    def check_turbine_location(self,x,y):
        new_x = np.array(x,dtype=float)
        new_y = np.array(y,dtype=float)
        moved = max(max(abs(new_x-self.old_x)),max(abs(new_y-self.old_y)))
        # print(moved,self.move_tol, moved > self.move_tol)
        if moved > self.move_tol:
            self.sparse_ids = self.construct_sparse_ids(self.x, self.sparse_RDs)
            self.vector_sparse_ids()
        self.old_x = new_x
        self.old_y = new_y

    def vector_sparse_ids(self):
        self.all_ids = np.zeros(len(self.sparse_ids)*self.dim,dtype=int)
        self.all_ids[0::self.dim] = self.dim*(self.sparse_ids)+0
        self.all_ids[1::self.dim] = self.dim*(self.sparse_ids)+1
        if self.dim == 3:
            self.all_ids[2::self.dim] = self.dim*(self.sparse_ids)+2
        self.all_ids = list(self.all_ids)

    def prepare_recompute_component(self, inputs, relevant_outputs):
        ### update the new controls inside the windfarm ###
        x = inputs[0]
        y = inputs[1]
        yaw = inputs[2]
        a = inputs[3]

        ### determine if we need to recalculate sparse_ids ###
        self.check_turbine_location(x,y)

        prepared = self.turb.build_actuator_disk(self.x,self.inflow_angle,self.sparse_ids)

        ### Update the turbine's actuator disk ###
        self.turb.actuator_disk = prepared

        return prepared

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        ### update the new controls inside the windfarm ###
        x = inputs[0]
        y = inputs[1]
        yaw = inputs[2]
        a = inputs[3]

        for test in adj_inputs:
            print(dir(test))
            print(dir(test.block_variable))
            print(test.block_variable.tag)
            exit()

        ### determine if we need to recalculate sparse_ids ###
        self.check_turbine_location(x,y)

        prepared = {}
        for name in self.control_types:
            d_actuator = self.turb.build_actuator_disk(self.x,self.inflow_angle,self.sparse_ids,dfd=name)
            prepared[name] = d_actuator.vector().get_local()

        return prepared

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name = block_variable.tag
        print("Calculating Derivative: " +name+"_"+repr(self.turb.index))

        ### Apply derivative to previous in tape ###
        adj_output = np.inner(adj_inputs[0].get_local(self.all_ids),prepared[name])

        ### MPI communication
        adj_output = np.float64(adj_output)
        recv_buff = np.zeros(1, dtype=np.float64)
        self.farm.params.comm.Allreduce(adj_output, recv_buff)
        adj_output = recv_buff

        return np.array(adj_output)
