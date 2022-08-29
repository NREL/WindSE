from pyadjoint.block import Block
from pyadjoint.tape import no_annotations
from . import assemble, dx
import copy

class ControlUpdaterBlock(Block):
    def __init__(self, J, problem, **kwargs):
        super(ControlUpdaterBlock, self).__init__()
        self.problem = problem
        self.farm = problem.farm
        self.time = kwargs.get("time",None)
        self.problem.df_first_save = True
        self.dt_sum = copy.copy(self.problem.dt_sum)

        # Add dependencies on the controls
        self.num_dependancies = 0
        for turb in self.farm.turbines:
            for control_name in turb.controls_list:
                control = getattr(turb,"m"+control_name)
                if isinstance(control,list):
                    for i in range(len(control)):
                        self.add_dependency(control[i])
                        control[i].block_variable.cu_tag = (control_name, turb.index, i)
                        self.num_dependancies += 1
                else:
                    self.add_dependency(control)
                    control.block_variable.cu_tag = (control_name, turb.index, -1)
                    self.num_dependancies += 1

        # self.add_dependency(self.problem.farm.turbines[0].simTime)
        # self.problem.farm.turbines[0].simTime.cu_tag = ("simTime",-1,-1)
        # self.num_dependancies += 1

        # self.add_dependency(self.problem.farm.turbines[0].simTime_prev)
        # self.problem.farm.turbines[0].simTime_prev.cu_tag = ("simTime_prev",-1,-1)
        # self.num_dependancies += 1


        ### TODO: we need to make this work with unsteady as well since unsteady doesn't use problem.up_k. possible solution just update u_k,p_k but check if len(u_k.vector()[:]) = len(inputs[i].vector()[:]) first
        self.add_dependency(self.problem.up_k)
        self.problem.up_k.block_variable.cu_tag = ("up",-1,-1)
        self.num_dependancies += 1

        self.add_dependency(self.problem.u_k)
        self.problem.u_k.block_variable.cu_tag = ("u",-1,-1)
        self.num_dependancies += 1

        self.add_dependency(self.problem.p_k)
        self.problem.p_k.block_variable.cu_tag = ("p",-1,-1)
        self.num_dependancies += 1

        J.block_variable.cu_tag = ("J",-1,-1)
        self.add_dependency(J)
        self.num_dependancies += 1

    def __str__(self):
        return "ControlUpdaterBlock"

    def recompute_component(self, inputs, block_variable, idx, prepared):
        # print(f"simTime =      {float(inputs[-6])}")
        # print(f"simTime_prev = {float(inputs[-5])}")
        print(f"Current Objective = {inputs[-1]/self.dt_sum if self.dt_sum>0 else inputs[-1]}")
        print(f"u(0, 0,150):   {inputs[-3]([0.0, 0.0,150.0])}")
        print(f"u(0, 0,210):   {inputs[-3]([0.0, 0.0,210.0])}")
        print(f"u(0,60,150):   {inputs[-3]([0.0,60.0,150.0])}")
        print(f"max(u):        {inputs[-3].vector().max()}")
        print(f"min(u):        {inputs[-3].vector().min()}")
        print(f"integral(u_x): {assemble(inputs[-3][0]*dx)}")
        # if self.problem.df_first_save:
        #     self.problem.df_velocity_file = self.problem.params.Save(inputs[-3],"df_velocity",subfolder="timeSeries/",val=self.time)
        #     self.problem.df_first_save = False
        # else:
        #     self.problem.params.Save(inputs[-3],"df_velocity",subfolder="timeSeries/",val=self.time,file=self.problem.df_velocity_file)

        return inputs[-1]

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        name, turb_id, seg_id = block_variable.cu_tag

        if name == "J":
            return adj_inputs[0]
        else:
            return 0.0

    def recompute(self, markings=False):
        print(f"Forward Solve Time: {self.time}")
        # self.Update()

        # Run original recompute
        Block.recompute(self, markings)
        return

    @no_annotations
    def evaluate_adj(self, markings=False):
        print(f"Adjoint Solve Time: {self.time}")
        # self.Update()

        # Run original evaluate_adj
        Block.evaluate_adj(self, markings)
        return

    def Update(self):
        # Get new dependency values
        deps = self.get_dependencies()
        inputs = [bv.saved_output for bv in deps]

        # update farm floats
        for i in range(len(inputs)):
            name, turb_id, seg_id = deps[i].cu_tag
            if turb_id != -1:
                if seg_id != -1:
                    control = getattr(self.farm.turbines[turb_id],name)
                    control[seg_id] = float(inputs[i])
                else:
                    setattr(self.farm.turbines[turb_id],name,float(inputs[i]))
            elif name in ["up", "u", "k"]:
                self.problem.up_k.assign(inputs[i])
        self.farm.update_controls()
