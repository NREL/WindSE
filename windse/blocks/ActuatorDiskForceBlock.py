from pyadjoint.block import Block

class ActuatorDiskForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x,farm,fs,sparse_ids,**kwargs):
        super(ActuatorDiskForceBlock, self).__init__()
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
            self.farm.mx[i].block_variable.tag = ("x",i,-1)
            self.add_dependency(self.farm.mx[i])

            self.farm.my[i].block_variable.tag = ("y",i,-1)
            self.add_dependency(self.farm.my[i])

            self.farm.myaw[i].block_variable.tag = ("yaw",i,-1)
            self.add_dependency(self.farm.myaw[i])

            self.farm.ma[i].block_variable.tag = ("a",i,-1)
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
        return "ActuatorLineForceBlock"

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
        name, turb_idx, _ = block_variable.tag
        print("Calculating Derivative: " +name+"_"+repr(turb_idx))
        ### Apply derivative to previous in tape ###
        adj_output = 0
        for i in range(3):
            adj_output += np.inner(adj_inputs[i].get_local(self.all_ids),prepared[name][i][:,turb_idx])

        adj_output = np.float64(adj_output)
        recv_buff = np.zeros(1, dtype=np.float64)
        self.farm.params.comm.Allreduce(adj_output, recv_buff)
        adj_output = recv_buff

        # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
        return np.array(adj_output)
