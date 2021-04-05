import dolfin
import dolfin_adjoint
import numpy as np
from sys import platform
import os,shutil
import time
import copy
import ufl

from windse.helper_functions import BaseHeight as backend_BaseHeight
from windse.helper_functions import RadialChordForce as backend_RadialChordForce
from windse.helper_functions import CalculateDiskTurbineForces as backend_CalculateDiskTurbineForces
from windse.helper_functions import UpdateActuatorLineForce as backend_UpdateActuatorLineForce
from windse import windse_parameters
from pyadjoint.tape import get_working_tape, annotate_tape, stop_annotating, no_annotations
from pyadjoint.block import Block
from pyadjoint.overloaded_type import create_overloaded_object
from pyadjoint.enlisting import Enlist
import cProfile

### This import improves the plotter functionality on Mac ###
if platform == 'darwin':
    import matplotlib
    matplotlib.use('TKAgg')
import matplotlib.pyplot as plt

def BaseHeight(x,y,ground,**kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)
    with stop_annotating():
        output = backend_BaseHeight(x,y,ground)
    output = create_overloaded_object(output)

    if annotate:        
        block = BaseHeightBlock(x, y, ground, output)
        block.add_output(output.create_block_variable())

        tape = get_working_tape()
        tape.add_block(block)

    return output

class BaseHeightBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, x, y, ground, output):
        super(BaseHeightBlock, self).__init__()
        self.x = x
        self.y = y
        self.ground = ground
        self.output = output
        self.add_dependency(x)
        self.add_dependency(y)

    def __str__(self):
        return "BaseHeightBlock"

    def prepare_recompute_component(self, inputs, relevant_outputs):
        x = inputs[0]
        y = inputs[1]
        return [x, y]

    def recompute_component(self, inputs, block_variable, idx, prepared):
        x = prepared[0]
        y = prepared[1]
        return backend_BaseHeight(x,y,self.ground)

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        x = inputs[0]
        y = inputs[1]
        return [x, y]

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        x = prepared[0]
        y = prepared[1]
        # h = dolfin.Constant(10)
        adj_input = adj_inputs[0]

        # print((float(x),float(y)))
        # print(self.ground(x,y))
        # print((self.ground(x,y,dx=1),self.ground(x,y,dy=1)))

        if idx == 0:
            # print("taking x derivative")
            adj = self.ground(x,y,dx=1)
            # adj = (self.ground(x+h,y)-self.ground(x-h,y))/(2*h)
        elif idx == 1:
            # print("taking y derivative")
            adj = self.ground(x,y,dy=1)
            # adj = (self.ground(x,y+h)-self.ground(x,y-h))/(2*h)


        return adj_input * adj


def RadialChordForce(r,chord):

    return backend_RadialChordForce(r,chord)


# def ReynoldsSum(u,idx):
    
#     return backend_FunctionAdd(x1)

# def






def CalculateDiskTurbineForces(x,farm,fs,**kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)

    ### Get the actual output ###
    with stop_annotating():
        func_output, sparse_ids, actuator_array = backend_CalculateDiskTurbineForces(x,farm,fs,**kwargs)

    ### Create the block object ###
    if annotate:
        block = ActuatorDiskForceBlock(x,farm,fs,sparse_ids,**kwargs)

    ### Prepare outputs for recording on tape ###
    r = []
    func_output = Enlist(func_output)
    for out in func_output:
        if not isinstance(out,list):
            r.append(create_overloaded_object(out))

    ### Add everything to tape ###
    if annotate:
        for out in r:
            if not isinstance(out,list):
                block.add_output(out.create_block_variable())
        tape = get_working_tape()
        tape.add_block(block)

    return func_output.delist(r), sparse_ids, actuator_array

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

        return np.array(adj_output)





def Marker(u, simTime, out_file, **kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)
    u = create_overloaded_object(u)

    if annotate:        
        block = MarkerBlock(u, simTime, out_file)
        block.add_output(u.create_block_variable())

        tape = get_working_tape()
        tape.add_block(block)

    return u

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


def ControlUpdater(J ,problem, **kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)
    J = create_overloaded_object(J)

    if annotate:        
        block = ControlUpdaterBlock(J, problem)
        block.add_output(J.create_block_variable())

        tape = get_working_tape()
        tape.add_block(block)

    return J

class ControlUpdaterBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, J, problem):
        super(ControlUpdaterBlock, self).__init__()
        self.problem = problem
        self.farm = problem.farm

        self.control_dict = {
                                "c_lift": self.farm.cl,
                                "c_drag": self.farm.cd,
                                "chord": self.farm.chord,
                                "x": self.farm.x,
                                "y": self.farm.y,
                                "yaw": self.farm.yaw,
                                "a": self.farm.axial,
        }

        # Add dependencies on the controls
        self.num_dependancies = 0
        for i in range(self.problem.farm.numturbs):
            if self.farm.turbine_method == "alm" or self.farm.force == "chord": 
                for j in range(self.problem.num_blade_segments):
                    self.farm.mcl[i][j].block_variable.tag = ("c_lift", i, j)
                    self.add_dependency(self.farm.mcl[i][j])
                    self.num_dependancies += 1

                    self.farm.mcd[i][j].block_variable.tag = ("c_drag", i, j)
                    self.add_dependency(self.farm.mcd[i][j])
                    self.num_dependancies += 1

                    self.farm.mchord[i][j].block_variable.tag = ("chord", i, j)
                    self.add_dependency(self.farm.mchord[i][j])
                    self.num_dependancies += 1

            self.farm.mx[i].block_variable.tag = ("x",i,-1)
            self.add_dependency(self.farm.mx[i])
            self.num_dependancies += 1

            self.farm.my[i].block_variable.tag = ("y",i,-1)
            self.add_dependency(self.farm.my[i])
            self.num_dependancies += 1

            self.farm.myaw[i].block_variable.tag = ("yaw",i,-1)
            self.add_dependency(self.farm.myaw[i])
            self.num_dependancies += 1

            self.farm.ma[i].block_variable.tag = ("a",i,-1)
            self.add_dependency(self.farm.ma[i])
            self.num_dependancies += 1

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
            self.farm.SimpleControlUpdate()
            return adj_inputs[0]
        elif name in ["c_lift","c_drag","chord"]:
            self.control_dict[name][turb_id][seg_id] = float(inputs[idx])
            return 0.0
        else:
            self.control_dict[name][turb_id] = float(inputs[idx])
            return 0.0















# ================================================================
# ================================================================

# def UpdateActuatorLineForce(problem, u_local, simTime_id, turb_i, dt, **kwargs):
def UpdateActuatorLineForce(problem, mpi_u_fluid_constant, simTime_id, dt, turb_i, **kwargs):
    '''This is the adjoint version of RelativeHeight. It's goal is to 
    calculate the height of the turbine's base. At the same time, it creates
    a block that helps propagate the adjoint information.'''
    annotate = annotate_tape(kwargs)

    # Example of how to call the update actuator force function
    # UpdateActuatorLineForce(problem, simTime, dfd, tf):

    ### Get the actual output ###
    with stop_annotating():
        func_output = backend_UpdateActuatorLineForce(problem, mpi_u_fluid_constant, simTime_id, dt, turb_i, **kwargs)


    ### Create the block object ###
    if annotate:
        block = ActuatorLineForceBlock(problem, mpi_u_fluid_constant, simTime_id, dt, turb_i, **kwargs)

    ### Prepare outputs for recording on tape ###
    r = []
    func_output = Enlist(func_output)
    for out in func_output:
        if not isinstance(out,list):
            r.append(create_overloaded_object(out))

    ### Add everything to tape ###
    if annotate:
        for out in r:
            if not isinstance(out,list):
                block.add_output(out.create_block_variable())
        tape = get_working_tape()
        tape.add_block(block)

    return func_output.delist(r)

# ================================================================
# ================================================================

class ActuatorLineForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, problem, mpi_u_fluid_constant, simTime_id, dt, turb_i, **kwargs):
        super(ActuatorLineForceBlock, self).__init__()
        # self.problem = copy.copy(problem)
        # self.simTime = copy.copy(simTime)
        # self.turb_i = copy.copy(turb_i)
        # self.u_local = u_local.copy()

        self.problem = problem
        self.simTime_id = copy.copy(simTime_id)
        self.simTime = self.problem.simTime_list[simTime_id]
        self.dt = dt
        self.turb_i = turb_i
        self.u_local = problem.u_k
        # self.mpi_u_fluid_constant = mpi_u_fluid_constant

        # Add dependencies on the controls
        for j in range(self.problem.num_blade_segments):
            self.problem.mcl[self.turb_i][j].block_variable.tag = ("c_lift", self.turb_i, j)
            self.add_dependency(self.problem.mcl[self.turb_i][j])

            self.problem.mcd[self.turb_i][j].block_variable.tag = ("c_drag", self.turb_i, j)
            self.add_dependency(self.problem.mcd[self.turb_i][j])

            self.problem.mchord[self.turb_i][j].block_variable.tag = ("chord", self.turb_i, j)
            self.add_dependency(self.problem.mchord[self.turb_i][j])


        self.problem.farm.myaw[self.turb_i].block_variable.tag = ("yaw", self.turb_i, -1)
        self.add_dependency(self.problem.farm.myaw[self.turb_i])

        self.problem.u_k1.block_variable.tag = ("u_local",self.turb_i,-1)
        self.add_dependency(self.problem.u_k1)

        # FIXME: This needs to be made a dolfin object so that it can be tracked
        self.problem.mpi_u_fluid_constant.block_variable.tag = ("mpi_u_fluid",self.turb_i,-1)
        self.add_dependency(self.problem.mpi_u_fluid_constant)

        # Tabulate which controls matter
        self.control_types = []
        if "lift" in problem.farm.control_types:
            self.control_types.append("c_lift")
        if "drag" in problem.farm.control_types:
            self.control_types.append("c_drag")
        if "chord" in problem.farm.control_types:
            self.control_types.append("chord")
        if "yaw" in problem.farm.control_types:
            self.control_types.append("yaw")

        # Record the current values for controls
        self.old_c_lift = np.array(self.problem.mcl, dtype=float)
        self.old_c_drag = np.array(self.problem.mcd, dtype=float)
        self.old_chord = np.array(self.problem.mchord, dtype=float)


    def __str__(self):
        return "ActuatorLineForceBlock"


    def prepare_recompute_component(self, inputs, relevant_outputs):
        # update the new controls inside the windfarm
        c_lift = inputs[0:-3:3]
        c_drag = inputs[1:-3:3]
        chord =  inputs[2:-3:3]
        yaw  = inputs[-3]
        u_k1 = inputs[-2]
        mpi_u_fluid = inputs[-1]

        # c_lift = inputs[0:-1:3]
        # c_drag = inputs[1:-1:3]
        # chord =  inputs[2:-1:3]
        # yaw  = inputs[-1]
        # u_k1 = self.u_local

        # print()
        # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        # print("current controls from dolfin_adjoint_helper")
        # print(np.array(chord,dtype=float))
        # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        # print()

        self.problem.fprint("Current Forward Time: "+repr(self.simTime),tab=1)
        self.problem.fprint("Turbine "+repr(self.turb_i),tab=2)
        self.problem.fprint("Current Yaw:   "+repr(float(self.problem.farm.myaw[self.turb_i])),tab=2)
        self.problem.fprint("Current Chord: "+str(np.array(self.problem.mchord[self.turb_i],dtype=float)),tab=2)
        # self.problem.fprint("",special="footer")


        self.problem.UpdateActuatorLineControls(c_lift = c_lift, c_drag = c_drag, chord = chord, yaw=yaw, turb_i=self.turb_i)

        # print()
        # print("uk("+repr(self.simTime)+")   = "+repr(np.mean(self.problem.u_k.vector()[:])))
        # print("uk1("+repr(self.simTime)+")  = "+repr(np.mean(self.problem.u_k1.vector()[:])))
        # print("uk2("+repr(self.simTime)+")  = "+repr(np.mean(self.problem.u_k2.vector()[:])))
        # print("uloc("+repr(self.simTime)+") = "+repr(np.mean(self.u_local.vector()[:])))
        # print("u_k1 = "+repr(np.mean(u_k1.vector()[:])))
        # print("uloc = "+repr(np.mean(self.u_local.vector()[:])))

        # Since dfd=None here, prepared is a dolfin function (tf) [1 x numPts*ndim]
        # prepared = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
        prepared = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid, self.simTime_id, self.dt, self.turb_i)
        # print("tf   = "+repr(np.mean(prepared.vector()[:])))
        # print()

        ### Calculate the exact current value of the objective function
        ### ONLY WORKS WITH ALM POWER OBJECTIVE AND ONE TURBINE 
        if self.simTime>=self.problem.record_time and self.problem.farm.numturbs==1:
            t_ids = np.logical_and(np.array(self.problem.simTime_list)>=self.problem.record_time,np.array(self.problem.simTime_list)<=self.simTime)
            dts = np.array(self.problem.dt_list)[t_ids]
            tos = np.array(self.problem.rotor_torque_dolfin_time)[t_ids]
            val = dts*tos 
            self.obj_value = (2.0*np.pi*self.problem.rpm/60.0)*np.sum(val)/np.sum(dts)/1.0e6
            self.problem.fprint("Current Power: "+repr(self.obj_value),tab=2)


        return prepared

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        ### update the new controls inside the windfarm ###
        c_lift = inputs[0:-3:3]
        c_drag = inputs[1:-3:3]
        chord =  inputs[2:-3:3]
        yaw  = inputs[-3]
        u_k1 = inputs[-2]
        mpi_u_fluid = inputs[-1]

        # c_lift = inputs[0:-1:3]
        # c_drag = inputs[1:-1:3]
        # chord =  inputs[2:-1:3]
        # yaw  = inputs[-1]
        # u_k1 = self.u_local

        self.problem.fprint("Current Adjoint Time: "+repr(self.simTime),tab=1)
        self.problem.fprint("Turbine "+repr(self.turb_i),tab=2)
        self.problem.fprint("Current Yaw:   "+repr(float(self.problem.farm.myaw[self.turb_i])),tab=2)
        self.problem.fprint("Current Chord: "+str(np.array(self.problem.mchord[self.turb_i],dtype=float)),tab=2)
        if self.simTime>=self.problem.record_time and self.problem.farm.numturbs==1:
            t_ids = np.logical_and(np.array(self.problem.simTime_list)>=self.problem.record_time,np.array(self.problem.simTime_list)<=self.simTime)
            dts = np.array(self.problem.dt_list)[t_ids]
            tos = np.array(self.problem.rotor_torque_dolfin_time)[t_ids]
            val = dts*tos 
            self.obj_value = (2.0*np.pi*self.problem.rpm/60.0)*np.sum(val)/np.sum(dts)/1.0e6
            self.problem.fprint("Current Power: "+repr(self.obj_value),tab=2)
        # self.problem.fprint("",special="footer")


        prepared = {}

        self.problem.UpdateActuatorLineControls(c_lift = c_lift, c_drag = c_drag, chord = chord, yaw=yaw, turb_i=self.turb_i)

        # for name in self.control_types:
        #     # pr = cProfile.Profile()
        #     # pr.enable()
        #     # print("calculating "+name+" derivatives")
        #     # Since dfd is not None here, prepared is a Numpy array of derivatives [numPts*ndim x numControls] 
        #     prepared[name] = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i, dfd=name, verbose=True)
        #     # pr.disable()
        #     # pr.print_stats()

        #     # exit()


        if "chord" in self.control_types:
            prepared["chord"] = []
            for i in range(self.problem.num_blade_segments):
                old_chord_value = copy.copy(chord[i])
                h_mag = 0.0001#*old_chord_value
                chord[i] = old_chord_value+h_mag
                self.problem.UpdateActuatorLineControls(chord = chord, turb_i=self.turb_i)
                # temp_uph = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
                temp_uph = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid, self.simTime_id, self.dt, self.turb_i)
                chord[i] = old_chord_value-h_mag
                self.problem.UpdateActuatorLineControls(chord = chord, turb_i=self.turb_i)
                # temp_umh = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
                temp_umh = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid, self.simTime_id, self.dt, self.turb_i)
                chord[i] = old_chord_value
                self.problem.UpdateActuatorLineControls(chord = chord, turb_i=self.turb_i)
                dtf_dc = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)
                prepared["chord"].append(dtf_dc)
            #FIXME: combine all these vectors onto rank 0, sum them, then broadcast that as the prepared vector
            prepared["chord"] = np.array(prepared["chord"]).T

        if "yaw" in self.control_types:
            prepared["yaw"] = []
            old_yaw_value = copy.copy(yaw)
            h_mag = 0.0001#*old_chord_value
            yaw = old_yaw_value+h_mag
            self.problem.UpdateActuatorLineControls(yaw = yaw, turb_i=self.turb_i)
            # temp_uph = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
            temp_uph = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid, self.simTime_id, self.dt, self.turb_i)
            yaw = old_yaw_value-h_mag
            self.problem.UpdateActuatorLineControls(yaw = yaw, turb_i=self.turb_i)
            # temp_umh = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
            temp_umh = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid, self.simTime_id, self.dt, self.turb_i)
            yaw = old_yaw_value
            self.problem.UpdateActuatorLineControls(yaw = yaw, turb_i=self.turb_i)
            dyaw_dc = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)
            # print(min(temp_uph.vector().get_local()),max(temp_uph.vector().get_local()))
            # print(min(temp_umh.vector().get_local()),max(temp_umh.vector().get_local()))
            prepared["yaw"] = dyaw_dc












        # tf = backend_UpdateActuatorLineForce(self.problem, u_k1, self.simTime_id, self.dt, self.turb_i)
        # gut = dolfin.grad(u_k1)
        # gtt = dolfin.grad(tf)
        # gu = gut#.T
        # gt = gtt#.T
        # prepared["u_local"] = [[],[],[]]



        # if not hasattr(self.problem,"dtfdu_files"):
        #     self.problem.dtfdu_files = [dolfin.File(self.problem.params.folder+"debug/dtf0.pvd"),
        #                                 dolfin.File(self.problem.params.folder+"debug/dtf1.pvd"),
        #                                 dolfin.File(self.problem.params.folder+"debug/dtf2.pvd")]


        # for i in range(3):

        #     test = dolfin.Function(self.problem.fs.V)
        #     test.rename("dtf"+repr(i),"dtf"+repr(i))
        #     test_vec = test.vector().get_local()

        #     # for j in range(3):
        #     #     this only work for 0 degree inflow (west to east or east to west)
        #     #     if self.simTime == 0:
        #     #         # val1 = gu[j,0]
        #     #         # val1 = dolfin.project(val1, self.problem.fs.Q,solver_type='cg', preconditioner_type='hypre_amg')

        #     #         # val2 = gu[j,1]
        #     #         # val2 = dolfin.project(val2, self.problem.fs.Q,solver_type='cg', preconditioner_type='hypre_amg')

        #     #         # val3 = gu[j,2]
        #     #         # val3 = dolfin.project(val3, self.problem.fs.Q,solver_type='cg', preconditioner_type='hypre_amg')

        #     #         # print(min(val1.vector().get_local()),min(val2.vector().get_local()),min(val3.vector().get_local()))
        #     #         # print(max(val1.vector().get_local()),max(val2.vector().get_local()),max(val3.vector().get_local()))
        #     #         # if j == 0:
        #     #         #     val = gt[i,0]*(1/gu[j,0])+gt[i,1]*(1/gu[j,1])+gt[i,2]*(1/gu[j,1])
        #     #         # else:
        #     #         val = dolfin.Constant(0.0)
        #     #     else:
        #     #         val = gt[i,0]*(1/gu[j,0])+gt[i,1]*(1/gu[j,1])+gt[i,2]*(1/gu[j,1])

        #     if self.simTime == 0:
        #         temp = dolfin.Function(self.problem.fs.Q)
        #         temp.vector()[:] = 0.0
        #         temp = temp.vector().get_local()
        #         for j in range(3):
        #             # val1 = gu[i,j]
        #             # val1 = dolfin.project(val1, self.problem.fs.Q,solver_type='cg', preconditioner_type='hypre_amg')
        #             # print("du_"+repr(i)+"/d_"+repr(j)+repr((min(val1.vector().get_local()),max(val1.vector().get_local()))))
        #             test_vec[j::3] = 0.0

        #             prepared["u_local"][i].append(temp)

        #     else:
        #         ### Create dTi/Dux ###
        #         val_numer =  gt[i,2]*gu[1,1]*gu[2,0] - gt[i,1]*gu[1,2]*gu[2,0] - gt[i,2]*gu[1,0]*gu[2,1] + gt[i,0]*gu[1,2]*gu[2,1] + gt[i,1]*gu[1,0]*gu[2,2] - gt[i,0]*gu[1,1]*gu[2,2]
        #         val_denom =  gu[0,2]*gu[1,1]*gu[2,0] - gu[0,1]*gu[1,2]*gu[2,0] - gu[0,2]*gu[1,0]*gu[2,1] + gu[0,0]*gu[1,2]*gu[2,1] + gu[0,1]*gu[1,0]*gu[2,2] - gu[0,0]*gu[1,1]*gu[2,2]
        #         val = val_numer/val_denom
        #         val = dolfin.project(val, self.problem.fs.Q,solver_type='mumps')
        #         test_vec[0::3] = val.vector().get_local()
        #         prepared["u_local"][i].append(val.vector().get_local())#/np.linalg.norm(val.vector().get_local()))

        #         ### Create dTi/Duy ###
        #         val_numer =  gt[i,2]*gu[0,1]*gu[2,0] - gt[i,1]*gu[0,2]*gu[2,0] - gt[i,2]*gu[0,0]*gu[2,1] + gt[i,0]*gu[0,2]*gu[2,1] + gt[i,1]*gu[0,0]*gu[2,2] - gt[i,0]*gu[0,1]*gu[2,2]
        #         val_denom = -gu[0,2]*gu[1,1]*gu[2,0] + gu[0,1]*gu[1,2]*gu[2,0] + gu[0,2]*gu[1,0]*gu[2,1] - gu[0,0]*gu[1,2]*gu[2,1] - gu[0,1]*gu[1,0]*gu[2,2] + gu[0,0]*gu[1,1]*gu[2,2]
        #         val = val_numer/val_denom
        #         val = dolfin.project(val, self.problem.fs.Q,solver_type='mumps')
        #         test_vec[1::3] = val.vector().get_local()
        #         prepared["u_local"][i].append(val.vector().get_local())#/np.linalg.norm(val.vector().get_local()))

        #         ### Create dTi/Duz ###
        #         val_numer =  gt[i,2]*gu[0,1]*gu[1,0] - gt[i,1]*gu[0,2]*gu[1,0] - gt[i,2]*gu[0,0]*gu[1,1] + gt[i,0]*gu[0,2]*gu[1,1] + gt[i,1]*gu[0,0]*gu[1,2] - gt[i,0]*gu[0,1]*gu[1,2]
        #         val_denom =  gu[0,2]*gu[1,1]*gu[2,0] - gu[0,1]*gu[1,2]*gu[2,0] - gu[0,2]*gu[1,0]*gu[2,1] + gu[0,0]*gu[1,2]*gu[2,1] + gu[0,1]*gu[1,0]*gu[2,2] - gu[0,0]*gu[1,1]*gu[2,2]
        #         val = val_numer/val_denom
        #         val = dolfin.project(val, self.problem.fs.Q,solver_type='mumps')
        #         test_vec[2::3] = val.vector().get_local()
        #         prepared["u_local"][i].append(val.vector().get_local())#/np.linalg.norm(val.vector().get_local()))

        #     test.vector()[:]=test_vec
        #     self.problem.dtfdu_files[i].write(test,self.simTime)





        if self.problem.farm.use_local_velocity:
            h_mag = 0.0001#*np.linalg.norm(u_local_vec)
            # u_local_vec = u_k1.vector().get_local()
            # print(np.linalg.norm(u_local_vec))

            mpi_u_fluid_buff = np.zeros(mpi_u_fluid.value_size())
            mpi_u_fluid.eval(mpi_u_fluid_buff, mpi_u_fluid_buff)

            prepared["u_local"] = []

            for i in range(self.problem.dom.dim):
                # h = np.zeros(u_local_vec.shape)
                # h[i::3] = h_mag
            
                # u_mh = dolfin.Function(u_k1.function_space())
                # u_mh.vector()[:] = u_local_vec-h


                # u_ph = dolfin.Function(u_k1.function_space())
                # u_ph.vector()[:] = u_local_vec+h

                # Perturb the mpi_u_fluid values in the x, y, and z-directions by +/- h
                h = np.zeros(mpi_u_fluid.value_size())
                h[i::3] = h_mag

                mpi_u_fluid_mh = dolfin.Constant(mpi_u_fluid_buff - h)
                mpi_u_fluid_ph = dolfin.Constant(mpi_u_fluid_buff + h)

                # temp_umh = backend_UpdateActuatorLineForce(self.problem, u_mh, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
                # temp_uph = backend_UpdateActuatorLineForce(self.problem, u_ph, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
                temp_umh = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid_mh, self.simTime_id, self.dt, self.turb_i)
                temp_uph = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid_ph, self.simTime_id, self.dt, self.turb_i)
                dtf_du = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)

                prepared["u_local"].append(dtf_du)

        return prepared

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, turbine_number, segment_index = block_variable.tag
        # print("calculating: " + name + "_" + repr(segment_index))


        # if not hasattr(self.problem,"ul_adjoint_files"):
        #     self.problem.ul_adjoint_files = dolfin.File(self.problem.params.folder+"debug/ul_adjoint.pvd")


        if name == "u_local":




            # adj_vec = adj_inputs[0].get_local()
            # comp_vec = np.zeros(adj_vec.shape)

            # for i in range(3):
            #     for j in range(3):
            #         comp_vec[i::3] += prepared[name][j][i]*adj_vec[j::3]

            # adj_output = dolfin.Function(self.u_local.function_space())
            # adj_output.vector()[:] = adj_vec #comp_vec
            # # adj_output.rename("u_adjoint","u_adjoint")
            # # self.problem.ul_adjoint_files.write(adj_output,self.simTime)
            





            # print(min(comp_vec),max(comp_vec))

            # print("I depend on u_local")
            #
            # prepared["u_local"] =  tf_x/u_x tf_y/u_x tf_z/u_x
            #                        tf_x/u_y tf_y/u_y tf_z/u_y
            #                        tf_x/u_z tf_y/u_z tf_z/u_z
            #
            # we need to calculate D = \nabla_u(tf)\cdot\nabla_c(u) which equals
            # 
            #    D_1 = tf_x/u_x * u_x/c + tf_x/u_x * u_x/c + tf_x/u_x * u_x/c
            #    D_2 = tf_y/u_x * u_x/c + tf_y/u_x * u_x/c + tf_y/u_x * u_x/c
            #    D_3 = tf_z/u_x * u_x/c + tf_z/u_x * u_x/c + tf_z/u_x * u_x/c
            #
            # u_j/c    = adj_vec[j::3]
            # tf_i/u_j = prepared["u_local"][j][i::3]
            # D_i      = comp_vec[i::3]
            #
            #
            #
            adj_output = dolfin.Function(self.u_local.function_space())
            if self.problem.farm.use_local_velocity:
                adj_vec = adj_inputs[0].get_local()
                comp_vec = np.zeros(adj_vec.shape)

                for i in range(3):
                    for j in range(3):
                        comp_vec[i::3] += prepared[name][j][i::3]*adj_vec[j::3]
                        # comp_vec[j::3] += prepared[name][j][i::3]*adj_vec[j::3]
                        # comp_vec[i::3] += np.inner(prepared[name][j][i::3],adj_vec[j::3])

                # for i in range(3):
                #     # comp_vec += prepared[name][i]*adj_vec
                #         comp_vec[i::3] += np.inner(prepared[name][j][i::3],adj_vec[j::3])

                adj_output.vector()[:] = comp_vec
                # adj_output.vector()[:] = 0.0
            else:
                adj_output.vector()[:] = 0.0


            # adj_output_vec = 0
            # for i in range(3):
            #     comp_vec = prepared[name][i]
            #     adj_vec = adj_inputs[0].get_local()
            #     adj_output_vec += np.inner(comp_vec, adj_vec)

            # adj_output = dolfin.Function(self.u_local.function_space())
            # adj_output.vector()[:] = adj_output_vec


            # adj_vec = adj_inputs[0].get_local()
            # adj_output = dolfin.Function(self.u_local.function_space())
            # adj_output.vector()[:] = (prepared[name][0]+prepared[name][1]+prepared[name][2])*adj_vec


            return adj_output.vector()
        elif name == "yaw":
            comp_vec = prepared[name]
            adj_vec = adj_inputs[0].get_local()
            adj_output = np.inner(comp_vec, adj_vec)
            return np.array(adj_output)
        else:
            comp_vec = prepared[name][:, segment_index]
            adj_vec = adj_inputs[0].get_local()
            adj_output = np.inner(comp_vec, adj_vec)
            return np.array(adj_output)


        ### Apply derivative to previous in tape ###

# ================================================================
# ================================================================
























#################################################################################################
############################## Overridden Dolfin_Adjoint Functions ##############################
#################################################################################################


### We need to override ALE move so that dolfin-adjoint can track the change ###
backend_move = dolfin.ALE.move
def move(mesh,bmesh, **kwargs):
    """ Refine is overloaded to ensure that the returned mesh is overloaded.
    """
    with stop_annotating():
        backend_move(mesh,bmesh, **kwargs)
    overloaded = create_overloaded_object(mesh)
    return overloaded

if "dolfin_adjoint_helper" not in dolfin.ALE.move.__module__:
    dolfin.ALE.move = move

def linalg_solve(*args, **kwargs):
    """This function overrides dolfin_adjoints.compat.linalg_solve.

    The original function doesn't allow for solver options because it uses
     the::

        dolfin.solve(A,x,b) 

    form which doesn't accept keyword arguments. However, It does except 
    additional arguments that defined some solver options, which we pass
    in manually

    Todo:

        Eventually, we want to replace this with a full PetscKrylovSolver()
        to get access to all the ksp options.

    """
    # print("performing a solve")
    return dolfin_adjoint.backend.solve(*args)#,"mumps") 

if "dolfin_adjoint_helper" not in dolfin_adjoint.types.compat.linalg_solve.__module__:
    dolfin_adjoint.types.compat.linalg_solve = linalg_solve

def assemble_adjoint_value(form, **kwargs):
    """Wrapper that assembles a matrix with boundary conditions"""
    bcs = kwargs.pop("bcs", ())
    # print(form)
    if windse_parameters["wind_farm"].get("turbine_method","dolfin") == "numpy":
        rep = 'tsfc'
    else:
        rep = 'uflacs'
    result = dolfin_adjoint.backend.assemble(form,form_compiler_parameters={'representation': rep})
    for bc in bcs:
        bc.apply(result)
    return result

if "dolfin_adjoint_helper" not in dolfin_adjoint.types.compat.assemble_adjoint_value.__module__:
    dolfin_adjoint.types.compat.assemble_adjoint_value = assemble_adjoint_value

shutil.rmtree(windse_parameters.folder+"debug/", ignore_errors=True)

def recompute_component(self, inputs, block_variable, idx, prepared):
    file_exists = False
    # print("resolving")

    # print()
    # print("hey Look at me I'm recomputing!")
    # print()
    """This function overrides 
    dolfin_adjoint.solving.SolveBlock.recompute_component

    The original function doesn't account for kwargs in a project so it 
    doesn't pass them to this function. For now, we just supply the
    solver_parameters manually in this case. 

    Todo:

        Eventually, we want to replace this with a full PetscKrylovSolver()
        to get access to all the ksp options.

    """
    lhs = prepared[0]
    rhs = prepared[1]
    func = prepared[2]
    bcs = prepared[3]


    if not self.forward_kwargs:
        dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, solver_parameters={'linear_solver': 'mumps'})
    else:
        dolfin_adjoint.backend.solve(lhs == rhs, func, bcs, **self.forward_kwargs)
        # print()
        # print("func = "+repr(np.mean(func.vector().get_local())))

        if "debug" in windse_parameters.output:
            if hasattr(self, 'solve_iteration'):
                self.solve_iteration += 1
            else:
                # print()
                # print("but we haven't been here before")
                self.recompute_set = 0
                while os.path.isfile(windse_parameters.folder+"debug/dolfin_adjoint_func_"+repr(self.recompute_set)+".pvd"):
                    self.recompute_set +=1
                # print("and that's the "+repr(self.recompute_set)+" time this has happened")
                # print()
                self.savefile = dolfin.File(windse_parameters.folder+"debug/dolfin_adjoint_func_"+repr(self.recompute_set)+".pvd")
                self.solve_iteration = 0
            
            
            u, p = func.split(True)
            u.rename("velocity","velocity")
            self.savefile << (u,self.solve_iteration)
    # print("assemble(func*dx): " + repr(float(dolfin.assemble(dolfin.inner(func,func)*dolfin.dx))))
    return func
if "dolfin_adjoint_helper" not in dolfin_adjoint.solving.SolveBlock.recompute_component.__module__:
    dolfin_adjoint.solving.SolveBlock.recompute_component = recompute_component

def _init_dependencies(self, *args, **kwargs):
    self.preconditioner_method = "default"
    self.solver_method = "mumps"

    if self.varform:
        eq = args[0]
        self.lhs = eq.lhs
        self.rhs = eq.rhs
        self.func = args[1]

        if len(args) > 2:
            self.bcs = args[2]
        elif "bcs" in kwargs:
            self.bcs = self.kwargs.pop("bcs")
            self.forward_kwargs.pop("bcs")
        else:
            self.bcs = []

        if self.bcs is None:
            self.bcs = []

        self.assemble_system = False
    else:
        # Linear algebra problem.
        # TODO: Consider checking if attributes exist.
        A = args[0]
        u = args[1]
        b = args[2]
        if len(args) >= 4:
            self.solver_method = args[3]
        if len(args) >= 5:
            self.preconditioner_method = args[4]

        sp = {"solver_parameters": {'linear_solver': self.solver_method,
                                    'preconditioner': self.preconditioner_method}}
        self.forward_kwargs.update(sp)

        self.lhs = A.form
        self.rhs = b.form
        self.bcs = A.bcs if hasattr(A, "bcs") else []
        self.func = u.function
        self.assemble_system = A.assemble_system if hasattr(A, "assemble_system") else False

    if not isinstance(self.bcs, list):
        self.bcs = [self.bcs]

    if isinstance(self.lhs, ufl.Form) and isinstance(self.rhs, ufl.Form):
        self.linear = True
        # Add dependence on coefficients on the right hand side.
        for c in self.rhs.coefficients():
            self.add_dependency(c, no_duplicates=True)
    else:
        self.linear = False

    for bc in self.bcs:
        self.add_dependency(bc, no_duplicates=True)

    for c in self.lhs.coefficients():
        self.add_dependency(c, no_duplicates=True)
if "dolfin_adjoint_helper" not in dolfin_adjoint.solving.SolveBlock._init_dependencies.__module__:
    dolfin_adjoint.solving.SolveBlock._init_dependencies = _init_dependencies

def _assemble_and_solve_adj_eq(self, dFdu_form, dJdu):
    dJdu_copy = dJdu.copy()
    kwargs = self.assemble_kwargs.copy()
    # Homogenize and apply boundary conditions on adj_dFdu and dJdu.
    bcs = self._homogenize_bcs()
    kwargs["bcs"] = bcs
    dFdu = dolfin_adjoint.compat.assemble_adjoint_value(dFdu_form, **kwargs)

    for bc in bcs:
        bc.apply(dJdu)

    adj_sol = dolfin_adjoint.Function(self.function_space)
    dolfin_adjoint.compat.linalg_solve(dFdu, adj_sol.vector(), dJdu, self.solver_method, self.preconditioner_method,**self.kwargs)

    adj_sol_bdy = dolfin_adjoint.compat.function_from_vector(self.function_space, dJdu_copy - dolfin_adjoint.compat.assemble_adjoint_value(
        dolfin.action(dFdu_form, adj_sol)))

    return adj_sol, adj_sol_bdy
if "dolfin_adjoint_helper" not in dolfin_adjoint.solving.SolveBlock._assemble_and_solve_adj_eq.__module__:
    dolfin_adjoint.solving.SolveBlock._assemble_and_solve_adj_eq = _assemble_and_solve_adj_eq





