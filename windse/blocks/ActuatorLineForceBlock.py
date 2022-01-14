import numpy as np
from pyadjoint.block import Block
import copy

from . import Function, Constant

class ActuatorLineForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, u_k, inflow_angle, fs, control_types=None, turb=None):
        super(ActuatorLineForceBlock, self).__init__()
        self.u_k = u_k
        self.inflow_angle = inflow_angle
        self.fs = fs
        self.control_types = control_types.copy()
        self.turb = turb

        self.DEBUGGING = False


        # Add dependencies on the controls
        for j in range(self.turb.num_blade_segments):
            self.turb.mcl[j].block_variable.tag = ("c_lift", self.turb.index, j)
            self.add_dependency(self.turb.mcl[j])

            self.turb.mcd[j].block_variable.tag = ("c_drag", self.turb.index, j)
            self.add_dependency(self.turb.mcd[j])

            self.turb.mchord[j].block_variable.tag = ("chord", self.turb.index, j)
            self.add_dependency(self.turb.mchord[j])


        self.turb.myaw.block_variable.tag = ("yaw", self.turb.index, -1)
        self.add_dependency(self.turb.myaw)

        self.turb.mpi_u_fluid_constant.block_variable.tag = ("mpi_u_fluid", self.turb.index, -1)
        self.add_dependency(self.turb.mpi_u_fluid_constant)

        if self.DEBUGGING:
            self.u_k.block_variable.tag = ("u_local",self.turb.index,-1)
            self.add_dependency(self.u_k)



    def __str__(self):
        return "ActuatorLineForceBlock"


    def prepare_recompute_component(self, inputs, relevant_outputs):
        if self.DEBUGGING:
            # update the new controls inside the windfarm
            self.turb.c_lift = np.array(inputs[0:-3:3], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-3:3], dtype=float)
            self.turb.chord =  np.array(inputs[2:-3:3], dtype=float)
            self.turb.yaw  = float(inputs[-3])
            self.turb.mpi_u_fluid = np.array(inputs[-2].values(), dtype=float)
            u_k = inputs[-1]

        else:
            self.turb.c_lift = np.array(inputs[0:-2:3], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-2:3], dtype=float)
            self.turb.chord =  np.array(inputs[2:-2:3], dtype=float)
            self.turb.yaw  = float(inputs[-2])
            self.turb.mpi_u_fluid = np.array(inputs[-1].values(), dtype=float)
            u_k = None

        self.turb.update_controls()

        self.turb.fprint("Current Forward Time: "+repr(self.turb.simTime),tab=1)
        self.turb.fprint("Turbine "+repr(self.turb.index),tab=2)
        self.turb.fprint("Current Yaw:   "+repr(float(self.turb.myaw)),tab=2)
        self.turb.fprint("Current Chord: "+str(np.array(self.turb.mchord,dtype=float)),tab=2)
        # self.turb.fprint("",special="footer")

        prepared = self.turb.build_actuator_lines(u_k, self.inflow_angle, self.fs)

        return prepared


    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared


    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        if self.DEBUGGING:
            # update the new controls inside the windfarm
            self.turb.c_lift = np.array(inputs[0:-3:3], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-3:3], dtype=float)
            self.turb.chord =  np.array(inputs[2:-3:3], dtype=float)
            self.turb.yaw  = float(inputs[-3])
            self.turb.mpi_u_fluid = np.array(inputs[-2].values(), dtype=float)
            u_k = inputs[-1]

        else:
            self.turb.c_lift = np.array(inputs[0:-2:3], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-2:3], dtype=float)
            self.turb.chord =  np.array(inputs[2:-2:3], dtype=float)
            self.turb.yaw  = float(inputs[-2])
            self.turb.mpi_u_fluid = np.array(inputs[-1].values(), dtype=float)
            u_k = None

        self.turb.update_controls()

        self.turb.fprint("Current Forward Time: "+repr(self.turb.simTime),tab=1)
        self.turb.fprint("Turbine "+repr(self.turb.index),tab=2)
        self.turb.fprint("Current Yaw:   "+repr(float(self.turb.myaw)),tab=2)
        self.turb.fprint("Current Chord: "+str(np.array(self.turb.mchord,dtype=float)),tab=2)

        prepared = {}

        if "chord" in self.control_types:
            prepared["chord"] = []
            for i in range(self.turb.num_blade_segments):
                old_chord_value = copy.copy(self.turb.chord[i])
                h_mag = 0.0001#*old_chord_value

                self.turb.chord[i] = old_chord_value+h_mag
                self.turb.update_controls()
                temp_uph = self.turb.build_actuator_lines(u_k, self.inflow_angle, self.fs)

                self.turb.chord[i] = old_chord_value-h_mag
                self.turb.update_controls()
                temp_umh = self.turb.build_actuator_lines(u_k, self.inflow_angle, self.fs)

                self.turb.chord[i] = old_chord_value
                self.turb.update_controls()

                dtf_dc = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)
                prepared["chord"].append(dtf_dc)

            #FIXME: combine all these vectors onto rank 0, sum them, then broadcast that as the prepared vector
            prepared["chord"] = np.array(prepared["chord"]).T


        if "yaw" in self.control_types:
            old_yaw_value = copy.copy(self.turb.yaw)
            h_mag = 0.0001#*old_chord_value

            self.turb.yaw = old_yaw_value+h_mag
            self.turb.update_controls()
            temp_uph = build_actuator_lines(u_k, self.inflow_angle, self.fs)

            self.turb.yaw = old_yaw_value-h_mag
            self.turb.update_controls()
            temp_umh = build_actuator_lines(u_k, self.inflow_angle, self.fs)

            self.turb.yaw = old_yaw_value
            self.turb.update_controls()

            dyaw_dc = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)
            prepared["yaw"] = dyaw_dc


        if self.turb.use_local_velocity:
            h_mag = 0.0001

            mpi_u_fluid_base = Constant(self.turb.mpi_u_fluid, name="mpi_u_fluid")

            prepared["u_local"] = []

            for i in range(self.turb.ndim):
                # Perturb the mpi_u_fluid values in the x, y, and z-directions by +/- h
                h = np.zeros(np.size(self.turb.mpi_u_fluid))
                h[i::self.turb.ndim] = h_mag

                mpi_u_fluid_ph = Constant(self.turb.mpi_u_fluid + h, name="mpi_u_fluid")
                mpi_u_fluid_mh = Constant(self.turb.mpi_u_fluid - h, name="mpi_u_fluid")

                self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_ph)
                temp_uph = self.turb.build_actuator_lines(u_k, self.inflow_angle, self.fs)

                self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_mh)
                temp_umh = self.turb.build_actuator_lines(u_k, self.inflow_angle, self.fs)

                self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_base)         
                dtf_du = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)

                prepared["u_local"].append(dtf_du)

        return prepared


    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, turbine_number, segment_index = block_variable.tag
        # print("calculating: " + name + "_" + repr(segment_index))
        
        print('evaluate_adj_component: %s' % (name))

        if name == "u_local":

            adj_output = Function(self.turb.tf.function_space())

            if self.turb.use_local_velocity:
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

            return adj_output.vector()

        elif name == "yaw":
            comp_vec = prepared[name]
            adj_vec = adj_inputs[0].get_local()
            adj_output = np.float64(np.inner(comp_vec, adj_vec))

            recv_buff = np.zeros(1, dtype=np.float64)
            self.turb.params.comm.Allreduce(adj_output, recv_buff)
            adj_output = recv_buff

            # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
            return np.array(adj_output)

        else:
            comp_vec = prepared[name][:, segment_index]
            adj_vec = adj_inputs[0].get_local()
            adj_output = np.float64(np.inner(comp_vec, adj_vec))

            recv_buff = np.zeros(1, dtype=np.float64)
            self.turb.params.comm.Allreduce(adj_output, recv_buff)
            adj_output = recv_buff

            # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
            return np.array(adj_output)
