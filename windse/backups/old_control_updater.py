OLD CONTROLS UPDATER BLOCK


def ControlUpdater(J ,problem, **kwargs):

    annotate = annotate_tape(kwargs)
    J = create_overloaded_object(J)

    if annotate:        
        block = ControlUpdaterBlock(J, problem,**kwargs)
        block.add_output(J.create_block_variable())

        tape = get_working_tape()
        tape.add_block(block)

    return J

class ControlUpdaterBlock(Block):
    def __init__(self, J, problem, **kwargs):
        super(ControlUpdaterBlock, self).__init__()
        self.problem = problem
        self.farm = problem.farm
        self.time = kwargs.get("time",None)
        self.control_dict = {
                                "c_lift": self.farm.cl,
                                "c_drag": self.farm.cd,
                                "chord": self.farm.chord,
                                "x": self.farm.x,
                                "y": self.farm.y,
                                "yaw": self.farm.yaw,
                                "a": self.farm.axial,
                                "up": self.problem.up_k,
        }

        # Add dependencies on the controls
        self.num_dependancies = 0
        for i in range(self.problem.farm.numturbs):
            if self.farm.turbine_method == "alm" or self.farm.force == "chord": 
                for j in range(self.problem.num_blade_segments):
                    self.add_dependency(self.farm.mcl[i][j])
                    self.farm.mcl[i][j].block_variable.tag = ("c_lift", i, j)
                    self.num_dependancies += 1

                    self.add_dependency(self.farm.mcd[i][j])
                    self.farm.mcd[i][j].block_variable.tag = ("c_drag", i, j)
                    self.num_dependancies += 1

                    self.add_dependency(self.farm.mchord[i][j])
                    self.farm.mchord[i][j].block_variable.tag = ("chord", i, j)
                    self.num_dependancies += 1

            self.add_dependency(self.farm.mx[i])
            self.farm.mx[i].block_variable.tag = ("x",i,-1)
            self.num_dependancies += 1

            self.add_dependency(self.farm.my[i])
            self.farm.my[i].block_variable.tag = ("y",i,-1)
            self.num_dependancies += 1

            self.add_dependency(self.farm.myaw[i])
            self.farm.myaw[i].block_variable.tag = ("yaw",i,-1)
            self.num_dependancies += 1

            self.add_dependency(self.farm.ma[i])
            self.farm.ma[i].block_variable.tag = ("a",i,-1)
            self.num_dependancies += 1

        ### TODO: we need to make this work with unsteady as well since unsteady doesn't use problem.up_k. possible solution just update u_k,p_k but check if len(u_k.vector()[:]) = len(inputs[i].vector()[:]) first
        self.add_dependency(self.problem.up_k)
        self.problem.up_k.block_variable.tag = ("up",-1,-1)
        self.num_dependancies += 1

        # self.add_dependency(self.problem.p_k)
        # self.problem.p_k.block_variable.tag = ("p",-1,-1)
        # self.num_dependancies += 1

        J.block_variable.tag = ("J",-1,-1)
        self.add_dependency(J)
        self.num_dependancies += 1

    def __str__(self):
        return "ControlUpdaterBlock"

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return inputs[-1]

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        name, turb_id, seg_id = block_variable.tag

        if name == "J":
            return adj_inputs[0]
        else:
            return 0.0

    def recompute(self, markings=False):
        print(f"Forward Solve Time: {self.time}")
        # with stop_annotating():
        self.Update()

        # Run original recompute
        Block.recompute(self, markings)
        return

    @no_annotations
    def evaluate_adj(self, markings=False):
        print(f"Adjoint Solve Time: {self.time}")
        # with stop_annotating():
        self.Update()

        # Run original evaluate_adj
        Block.evaluate_adj(self, markings)
        return

    def Update(self):
        # Get new dependency values
        deps = self.get_dependencies()
        inputs = [bv.saved_output for bv in deps]

        # update farm floats
        for i in range(len(inputs)):
            name, turb_id, seg_id = deps[i].tag
            if name != "J":
                if name in ["c_lift","c_drag","chord"]:
                    self.control_dict[name][turb_id][seg_id] = float(inputs[i])
                elif name in ["up"]:
                    # self.control_dict[name] = inputs[i]
                    self.control_dict[name].assign(inputs[i])
                    # self.control_dict[name].assign(inputs[i], annotate=False)
                else:
                    self.control_dict[name][turb_id] = float(inputs[i])

        # update farm Constants()
        # with stop_annotating():
        self.problem.SimpleControlUpdate()







    def CreateControls(self):

        ### Control pointers is a way of reference which parts of the original farm lists correspond to the controls. it is rather inefficient but it works so...
        self.controls = []
        self.control_pointers = []
        self.names = []
        self.indexes = [[],[],[],[],[],[],[]]
        self.init_vals = []
        j = 0

        if "layout" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[0].append(j)
                j+=1
                self.names.append("x_"+repr(i))
                self.controls.append(Control(self.farm.mx[i]))
                self.control_pointers.append((self.farm.x,i))
                self.init_vals.append(Constant(float(self.farm.mx[i])))

                self.indexes[1].append(j)
                j+=1
                self.names.append("y_"+repr(i))
                self.controls.append(Control(self.farm.my[i]))
                self.control_pointers.append((self.farm.y,i))
                self.init_vals.append(Constant(float(self.farm.my[i])))

        if "yaw" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[2].append(j)
                j+=1
                self.names.append("yaw_"+repr(i))
                self.controls.append(Control(self.farm.myaw[i]))
                self.control_pointers.append((self.farm.yaw,i))
                self.init_vals.append(Constant(float(self.farm.myaw[i])))

        if "axial" in self.control_types:
            for i in self.solver.opt_turb_id:
                self.indexes[3].append(j)
                j+=1
                self.names.append("axial_"+repr(i))
                self.controls.append(Control(self.farm.ma[i]))
                self.control_pointers.append((self.farm.a,i))
                self.init_vals.append(Constant(float(self.farm.ma[i])))

        if "lift" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    self.control_pointers.append((self.farm.cl,[i,k]))
                    self.indexes[4].append(j)
                    j+=1
                    self.names.append("lift_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.mcl[i][k]))
                    self.init_vals.append(Constant(float(self.farm.mcl[i][k])))

        if "drag" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    self.control_pointers.append((self.farm.cd,[i,k]))
                    self.indexes[5].append(j)
                    j+=1
                    self.names.append("drag_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.mcd[i][k]))
                    self.init_vals.append(Constant(float(self.farm.mcd[i][k])))

        if "chord" in self.control_types:
            for i in self.solver.opt_turb_id:
                for k in range(self.farm.num_blade_segments):
                    self.control_pointers.append((self.farm.chord,[i,k]))
                    self.indexes[6].append(j)
                    j+=1
                    self.names.append("chord_"+repr(i)+"_"+repr(k))
                    self.controls.append(Control(self.farm.mchord[i][k]))
                    self.init_vals.append(Constant(float(self.farm.mchord[i][k])))
        self.num_controls = len(self.controls)