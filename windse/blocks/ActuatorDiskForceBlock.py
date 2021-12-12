from pyadjoint.block import Block
import numpy as np

class ActuatorDiskForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x, inflow_angle, fs, sparse_ids, turb=None, construct_sparse_ids=None, actuator_disk=None, control_types=None):
        super(ActuatorDiskForceBlock, self).__init__()
        self.x = x
        self.fs = fs
        self.inflow_angle = inflow_angle
        self.sparse_ids = sparse_ids
        self.turb = turb
        self.construct_sparse_ids = construct_sparse_ids
        self.sparse_RDs = 2.0
        self.dim = turb.dom.dim
        self.control_types = control_types.copy()

        if "layout" in self.control_types:
            self.control_types += ["x","y"]
            self.control_types.remove("layout")
            
        ### Calculate how far the turbines need to move to recalculate sparse ids ###
        self.move_tol = 0.1*turb.dom.mesh.rmin()

        ### Add dependencies on the controls ###
        self.turb.mx.block_variable.tag = ("x",self.turb.index,-1)
        self.add_dependency(self.turb.mx)

        self.turb.my.block_variable.tag = ("y",self.turb.index,-1)
        self.add_dependency(self.turb.my)

        self.turb.myaw.block_variable.tag = ("yaw",self.turb.index,-1)
        self.add_dependency(self.turb.myaw)

        self.turb.maxial.block_variable.tag = ("axial",self.turb.index,-1)
        self.add_dependency(self.turb.maxial)

        ### Record the current values for x and y ###
        self.old_x = self.turb.x
        self.old_y = self.turb.x

        ### map ids to vector ids ###
        self.vector_sparse_ids()

    def __str__(self):
        return "ActuatorLineForceBlock"

    def check_turbine_location(self,x,y):
        new_x = float(x)
        new_y = float(y)
        # moved = np.sqrt(abs(new_x-self.old_x)**2 + abs(new_y-self.old_y)**2)
        moved = max(abs(new_x-self.old_x), abs(new_y-self.old_y))
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
        self.turb.x =     float(inputs[0])
        self.turb.y =     float(inputs[1])
        self.turb.yaw =   float(inputs[2])
        self.turb.axial = float(inputs[3])
        self.turb.update_controls()

        ### determine if we need to recalculate sparse_ids ###
        self.check_turbine_location(self.turb.x,self.turb.y)

        prepared = self.turb.build_actuator_disk(self.x,self.inflow_angle,self.fs,self.sparse_ids, actuator_disk=self.turb.actuator_disk)

        return prepared

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        print("Calculating derivatives for turbine: "+repr(self.turb.index))
        
        ### update the new controls inside the windfarm ###
        self.turb.x =     float(inputs[0])
        self.turb.y =     float(inputs[1])
        self.turb.yaw =   float(inputs[2])
        self.turb.axial = float(inputs[3])
        self.turb.update_controls()

        ### determine if we need to recalculate sparse_ids ###
        self.check_turbine_location(self.turb.x,self.turb.y)

        prepared = {}
        for name in self.control_types:
            prepared[name] = self.turb.build_actuator_disk(self.x,self.inflow_angle,self.fs,self.sparse_ids,dfd=name)

        return prepared

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, turb_id, _ = block_variable.tag

        ### Apply derivative to previous in tape ###
        adj_output = 0
        for i in range(3):
            adj_output += np.inner(adj_inputs[i].get_local(self.all_ids),prepared[name][i].T)

        ### MPI communication
        adj_output = np.float64(adj_output)
        recv_buff = np.zeros(1, dtype=np.float64)
        self.turb.params.comm.Allreduce(adj_output, recv_buff)
        adj_output = recv_buff

        return np.array(adj_output)
