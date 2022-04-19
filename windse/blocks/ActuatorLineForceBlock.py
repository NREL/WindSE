import numpy as np
from pyadjoint.block import Block
import copy

from . import Function, Constant

class ActuatorLineForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, u_k, inflow_angle, fs, control_types=None, turb=None, **kwargs):
        super(ActuatorLineForceBlock, self).__init__()
        self.u_k = u_k
        self.inflow_angle = inflow_angle
        self.fs = fs
        self.control_types = control_types.copy()
        self.turb = turb

        self.INCLUDE_U_K = True

        self.simTime = kwargs['simTime']
        self.simTime_prev = kwargs['simTime_prev']
        self.dt = kwargs['dt']

        # Add dependencies on the controls
        for j in range(self.turb.num_blade_segments):
            self.turb.mcl[j].block_variable.tag = ("c_lift", self.turb.index, j, -1)
            self.add_dependency(self.turb.mcl[j])

            self.turb.mcd[j].block_variable.tag = ("c_drag", self.turb.index, j, -1)
            self.add_dependency(self.turb.mcd[j])

            self.turb.mchord[j].block_variable.tag = ("chord", self.turb.index, j, -1)
            self.add_dependency(self.turb.mchord[j])


            for i in range(3):
                self.turb.mpi_u_fluid[i][j].block_variable.tag = ("mpi_u_fluid", self.turb.index, j, i)
                self.add_dependency(self.turb.mpi_u_fluid[i][j])

        self.turb.myaw.block_variable.tag = ("yaw", self.turb.index, -1, -1)
        self.add_dependency(self.turb.myaw)


        if self.INCLUDE_U_K:
            self.u_k.block_variable.tag = ("u_local",self.turb.index,-1, -1)
            self.add_dependency(self.u_k)



    def __str__(self):
        return "ActuatorLineForceBlock"


    def prepare_recompute_component(self, inputs, relevant_outputs):
        if self.INCLUDE_U_K:
            # update the new controls inside the windfarm
            self.turb.c_lift = np.array(inputs[0:-2:6], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-2:6], dtype=float)
            self.turb.chord =  np.array(inputs[2:-2:6], dtype=float)
            self.turb.mpi_u_fluid[0] = inputs[3:-2:6]
            self.turb.mpi_u_fluid[1] = inputs[4:-2:6]
            self.turb.mpi_u_fluid[2] = inputs[5:-2:6]
            self.turb.yaw  = float(inputs[-2])
            u_k = inputs[-1]

        else:
            self.turb.c_lift = np.array(inputs[0:-1:6], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-1:6], dtype=float)
            self.turb.chord =  np.array(inputs[2:-1:6], dtype=float)
            self.turb.mpi_u_fluid[0] = inputs[3:-2:6]
            self.turb.mpi_u_fluid[1] = inputs[4:-2:6]
            self.turb.mpi_u_fluid[2] = inputs[5:-2:6]
            self.turb.yaw  = float(inputs[-1])
            u_k = None

        self.turb.update_controls()

        self.turb.fprint("Current Forward Time: "+repr(self.turb.simTime),tab=1)
        self.turb.fprint("Turbine "+repr(self.turb.index),tab=2)
        self.turb.fprint("Current Yaw:   "+repr(float(self.turb.myaw)),tab=2)
        self.turb.fprint("Current Chord: "+str(np.array(self.turb.mchord,dtype=float)),tab=2)
        # self.turb.fprint("",special="footer")

        prepared = self.turb.turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

        return prepared


    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared


    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        if self.INCLUDE_U_K:
            # update the new controls inside the windfarm
            self.turb.c_lift = np.array(inputs[0:-2:6], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-2:6], dtype=float)
            self.turb.chord =  np.array(inputs[2:-2:6], dtype=float)
            self.turb.mpi_u_fluid[0] = inputs[3:-2:6]
            self.turb.mpi_u_fluid[1] = inputs[4:-2:6]
            self.turb.mpi_u_fluid[2] = inputs[5:-2:6]
            self.turb.yaw  = float(inputs[-2])
            u_k = inputs[-1]

        else:
            self.turb.c_lift = np.array(inputs[0:-1:6], dtype=float)
            self.turb.c_drag = np.array(inputs[1:-1:6], dtype=float)
            self.turb.chord =  np.array(inputs[2:-1:6], dtype=float)
            self.turb.mpi_u_fluid[0] = inputs[3:-2:6]
            self.turb.mpi_u_fluid[1] = inputs[4:-2:6]
            self.turb.mpi_u_fluid[2] = inputs[5:-2:6]
            self.turb.yaw  = float(inputs[-1])
            u_k = None

        self.turb.update_controls()

        self.turb.fprint("Current Adjoint Time: "+repr(self.turb.simTime),tab=1)
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
                temp_uph = self.turb.turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

                self.turb.chord[i] = old_chord_value-h_mag
                self.turb.update_controls()
                temp_umh = self.turb.turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

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
            temp_uph = turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

            self.turb.yaw = old_yaw_value-h_mag
            self.turb.update_controls()
            temp_umh = turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

            self.turb.yaw = old_yaw_value
            self.turb.update_controls()

            dyaw_dc = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)
            prepared["yaw"] = dyaw_dc


        # if self.turb.use_local_velocity:
        #     h_mag = 0.0001

        #     mpi_u_fluid_base = Constant(self.turb.mpi_u_fluid, name="mpi_u_fluid")

        #     prepared["u_local"] = []

        #     for i in range(self.turb.ndim):
        #         # Perturb the mpi_u_fluid values in the x, y, and z-directions by +/- h
        #         h = np.zeros(np.size(self.turb.mpi_u_fluid))
        #         h[i::self.turb.ndim] = h_mag

        #         mpi_u_fluid_ph = Constant(self.turb.mpi_u_fluid + h, name="mpi_u_fluid")
        #         mpi_u_fluid_mh = Constant(self.turb.mpi_u_fluid - h, name="mpi_u_fluid")

        #         self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_ph)
        #         temp_uph = self.turb.turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

        #         self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_mh)
        #         temp_umh = self.turb.turbine_force(u_k, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)

        #         self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_base)         
        #         dtf_du = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)

        #         # dtf_du[:] = 0.0
        #         prepared["u_local"].append(dtf_du)


        if self.turb.use_local_velocity:
            h_mag = 0.0001#*np.linalg.norm(u_local_vec)
            u_local_vec = u_k.vector().get_local()
            # print(np.linalg.norm(u_local_vec))

            # mpi_u_fluid_buff = np.zeros(mpi_u_fluid.value_size())
            # mpi_u_fluid.eval(mpi_u_fluid_buff, mpi_u_fluid_buff)

            prepared["u_local"] = []

            for i in range(self.turb.ndim):
                h = np.zeros(u_local_vec.shape)
                h[i::3] = h_mag
            
                u_mh = Function(u_k.function_space())
                u_mh.vector()[:] = u_local_vec-h


                u_ph = Function(u_k.function_space())
                u_ph.vector()[:] = u_local_vec+h

                # # Perturb the mpi_u_fluid values in the x, y, and z-directions by +/- h
                # h = np.zeros(mpi_u_fluid.value_size())
                # h[i::3] = h_mag

                # mpi_u_fluid_mh = dolfin.Constant(mpi_u_fluid_buff - h)
                # mpi_u_fluid_ph = dolfin.Constant(mpi_u_fluid_buff + h)

                # temp_umh = backend_UpdateActuatorLineForce(self.problem, u_mh, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
                # temp_uph = backend_UpdateActuatorLineForce(self.problem, u_ph, self.simTime_id, self.dt, self.turb_i, self.mpi_u_fluid)
                temp_uph = self.turb.turbine_force(u_ph, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)
                temp_umh = self.turb.turbine_force(u_mh, self.inflow_angle, self.fs,simTime=self.simTime,simTime_prev=self.simTime_prev,dt=self.dt)
                # temp_umh = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid_mh, self.simTime_id, self.dt, self.turb_i)
                # temp_uph = backend_UpdateActuatorLineForce(self.problem, mpi_u_fluid_ph, self.simTime_id, self.dt, self.turb_i)
                dtf_du = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)

                prepared["u_local"].append(dtf_du)

        return prepared


    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, turbine_number, segment_index, blade_id = block_variable.tag
        # print("calculating: " + name + "_" + repr(segment_index))
        # print("input_val: "+repr(adj_inputs[0].norm('l2')))

        # print('evaluate_adj_component: %s' % (name))

        if name == "u_local":

            adj_output = Function(self.turb.tf.function_space())

            if self.turb.use_local_velocity:
                adj_vec = adj_inputs[0].get_local()
                comp_vec = np.zeros(adj_vec.shape)

                for i in range(3):
                    for j in range(3):
                        # print('partial_val_norm (%d, %d): %e' % (i, j, np.linalg.norm(prepared[name][j][i::3])))
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

            # print("output_val: "+repr(adj_output.vector().norm('l2')))
            return adj_output.vector()

        elif name == "yaw":
            # print("input_val: "+repr(np.array(adj_output)))
            comp_vec = prepared[name]
            # print('yaw_comp_vec_norm: %e' % (np.linalg.norm(comp_vec)))

            adj_vec = adj_inputs[0].get_local()
            adj_output = np.float64(np.inner(comp_vec, adj_vec))

            recv_buff = np.zeros(1, dtype=np.float64)
            self.turb.params.comm.Allreduce(adj_output, recv_buff)
            adj_output = recv_buff

            # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
            # print("output_val: "+repr(float(adj_output)))
            return np.array(adj_output)

        else:
            comp_vec = prepared[name][:, segment_index]
            # print('else_comp_vec_norm: %e' % (np.linalg.norm(comp_vec)))
            adj_vec = adj_inputs[0].get_local()
            adj_output = np.float64(np.inner(comp_vec, adj_vec))

            recv_buff = np.zeros(1, dtype=np.float64)
            self.turb.params.comm.Allreduce(adj_output, recv_buff)
            adj_output = recv_buff

            # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
            # print("output_val: "+repr(float(adj_output)))
            return np.array(adj_output)
