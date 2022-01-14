import numpy as np
from pyadjoint.block import Block
import copy

from . import Function, Constant

class ActuatorLineForceBlock(Block):
    '''This is the Block class that will be used to calculate adjoint 
    information for optimizations. '''
    def __init__(self, fs, u_k, inflow_angle, control_types=None, turb=None):
        super(ActuatorLineForceBlock, self).__init__()
        # self.turb = copy.copy(problem)
        # self.simTime = copy.copy(simTime)
        # self.turb.index = copy.copy(turb_i)
        # self.u_local = u_local.copy()

        # self.turb = problem
        # self.simTime_id = copy.copy(simTime_id)
        # self.simTime = self.turb.simTime_list[simTime_id]
        # self.dt = dt
        # self.turb.index = turb_i
        # self.u_local = problem.u_k
        # self.mpi_u_fluid_constant = mpi_u_fluid_constant\
        self.fs = fs
        self.u_k = u_k
        self.inflow_angle = inflow_angle
        self.control_types = control_types.copy()
        self.turb = turb


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

        self.u_k.block_variable.tag = ("u_local",self.turb.index,-1)
        self.add_dependency(self.u_k)

        # FIXME: This needs to be made a dolfin object so that it can be tracked
        self.turb.mpi_u_fluid_constant.block_variable.tag = ("mpi_u_fluid", self.turb.index, -1)
        self.add_dependency(self.turb.mpi_u_fluid_constant)

        # # Tabulate which controls matter
        # self.control_types = []
        # if "lift" in problem.farm.control_types:
        #     self.control_types.append("c_lift")
        # if "drag" in problem.farm.control_types:
        #     self.control_types.append("c_drag")
        # if "chord" in problem.farm.control_types:
        #     self.control_types.append("chord")
        # if "yaw" in problem.farm.control_types:
        #     self.control_types.append("yaw")

        # Record the current values for controls
        self.old_c_lift = np.array(self.turb.mcl, dtype=float)
        self.old_c_drag = np.array(self.turb.mcd, dtype=float)
        self.old_chord = np.array(self.turb.mchord, dtype=float)
        self.old_twist= np.array(self.turb.mtwist, dtype=float)


    def __str__(self):
        return "ActuatorLineForceBlock"


    def prepare_recompute_component(self, inputs, relevant_outputs):
        # update the new controls inside the windfarm
        self.turb.c_lift = np.array(inputs[0:-3:3], dtype=float)
        self.turb.c_drag = np.array(inputs[1:-3:3], dtype=float)
        self.turb.chord =  np.array(inputs[2:-3:3], dtype=float)
        self.turb.yaw  = float(inputs[-3])
        u_k = inputs[-2]
        self.turb.mpi_u_fluid = inputs[-1]

        self.turb.update_controls()

        self.turb.fprint("Current Forward Time: "+repr(self.turb.simTime_prev),tab=1)
        self.turb.fprint("Turbine "+repr(self.turb.index),tab=2)
        self.turb.fprint("Current Yaw:   "+repr(float(self.turb.myaw)),tab=2)
        self.turb.fprint("Current Chord: "+str(np.array(self.turb.mchord,dtype=float)),tab=2)
        # self.turb.fprint("",special="footer")


        # self.turb.UpdateActuatorLineControls(c_lift = c_lift, c_drag = c_drag, chord = chord, yaw=yaw, turb_i=self.turb.index)


        prepared = self.turb.build_actuator_lines(self.fs, u_k, self.inflow_angle)


        ### Calculate the exact current value of the objective function
        ### ONLY WORKS WITH ALM POWER OBJECTIVE AND ONE TURBINE 
        # if self.simTime>=self.turb.record_time and self.turb.farm.numturbs==1:
        #     t_ids = np.logical_and(np.array(self.turb.simTime_list)>=self.turb.record_time,np.array(self.turb.simTime_list)<=self.simTime)
        #     dts = np.array(self.turb.dt_list)[t_ids]
        #     tos = np.array(self.turb.rotor_torque_dolfin_time)[t_ids]
        #     val = dts*tos 
        #     self.obj_value = (2.0*np.pi*self.turb.rpm/60.0)*np.sum(val)/np.sum(dts)/1.0e6
        #     self.turb.fprint("Current Power: "+repr(self.obj_value),tab=2)

        return prepared

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return prepared

    def prepare_evaluate_adj(self, inputs, adj_inputs, relevant_dependencies):
        ### update the new controls inside the windfarm ###
        self.turb.c_lift = np.array(inputs[0:-3:3], dtype=float)
        self.turb.c_drag = np.array(inputs[1:-3:3], dtype=float)
        self.turb.chord =  np.array(inputs[2:-3:3], dtype=float)
        self.turb.yaw  = float(inputs[-3])
        u_k = inputs[-2]
        self.turb.mpi_u_fluid = inputs[-1]

        self.turb.update_controls()

        self.turb.fprint("Current Forward Time: "+repr(self.turb.simTime_prev),tab=1)
        self.turb.fprint("Turbine "+repr(self.turb.index),tab=2)
        self.turb.fprint("Current Yaw:   "+repr(float(self.turb.myaw)),tab=2)
        self.turb.fprint("Current Chord: "+str(np.array(self.turb.mchord,dtype=float)),tab=2)


        # if self.simTime>=self.turb.record_time and self.turb.farm.numturbs==1:
        #     t_ids = np.logical_and(np.array(self.turb.simTime_list)>=self.turb.record_time,np.array(self.turb.simTime_list)<=self.simTime)
        #     dts = np.array(self.turb.dt_list)[t_ids]
        #     tos = np.array(self.turb.rotor_torque_dolfin_time)[t_ids]
        #     val = dts*tos 
        #     self.obj_value = (2.0*np.pi*self.turb.rpm/60.0)*np.sum(val)/np.sum(dts)/1.0e6
        #     self.turb.fprint("Current Power: "+repr(self.obj_value),tab=2)
        # self.turb.fprint("",special="footer")


        prepared = {}

        # self.turb.UpdateActuatorLineControls(c_lift = c_lift, c_drag = c_drag, chord = chord, yaw=yaw, turb_i=self.turb.index)

        # for name in self.control_types:
        #     # pr = cProfile.Profile()
        #     # pr.enable()
        #     # print("calculating "+name+" derivatives")
        #     # Since dfd is not None here, prepared is a Numpy array of derivatives [numPts*ndim x numControls] 
        #     prepared[name] = backend_UpdateActuatorLineForce(self.turb, u_k1, self.simTime_id, self.dt, self.turb.index, dfd=name, verbose=True)
        #     # pr.disable()
        #     # pr.print_stats()

        #     # exit()


        if "chord" in self.control_types:
            prepared["chord"] = []
            for i in range(self.turb.num_blade_segments):
                old_chord_value = copy.copy(self.turb.chord[i])
                h_mag = 0.0001#*old_chord_value

                self.turb.chord[i] = old_chord_value+h_mag
                self.turb.update_controls()
                temp_uph = self.turb.build_actuator_lines(self.fs, u_k, self.inflow_angle)

                self.turb.chord[i] = old_chord_value-h_mag
                self.turb.update_controls()
                temp_umh = self.turb.build_actuator_lines(self.fs, u_k, self.inflow_angle)

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
            temp_uph = build_actuator_lines(self.fs, u_k, self.inflow_angle)

            self.turb.yaw = old_yaw_value-h_mag
            self.turb.update_controls()
            temp_umh = build_actuator_lines(self.fs, u_k, self.inflow_angle)

            self.turb.yaw = old_yaw_value
            self.turb.update_controls()

            dyaw_dc = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)
            prepared["yaw"] = dyaw_dc


        if self.turb.use_local_velocity:
            h_mag = 0.0001#*np.linalg.norm(u_local_vec)
            # u_local_vec = u_k1.vector().get_local()
            # print(np.linalg.norm(u_local_vec))

            mpi_u_fluid_buff = np.zeros(self.turb.mpi_u_fluid_constant.value_size())
            self.turb.mpi_u_fluid.eval(mpi_u_fluid_buff, mpi_u_fluid_buff)

            prepared["u_local"] = []

            for i in range(self.turb.ndim):
                # h = np.zeros(u_local_vec.shape)
                # h[i::3] = h_mag
            
                # u_mh = dolfin.Function(u_k1.function_space())
                # u_mh.vector()[:] = u_local_vec-h


                # u_ph = dolfin.Function(u_k1.function_space())
                # u_ph.vector()[:] = u_local_vec+h

                # Perturb the mpi_u_fluid values in the x, y, and z-directions by +/- h
                h = np.zeros(self.turb.mpi_u_fluid.value_size())
                h[i::self.turb.ndim] = h_mag

                mpi_u_fluid_ph = Constant(mpi_u_fluid_buff + h)
                mpi_u_fluid_mh = Constant(mpi_u_fluid_buff - h)

                # self.turb.mpi_u_fluid_constant.assi = mpi_u_fluid_buff + h
                # self.turb.update_controls()


                # temp_umh = backend_UpdateActuatorLineForce(self.turb, u_mh, self.simTime_id, self.dt, self.turb.index, self.mpi_u_fluid)
                # temp_uph = backend_UpdateActuatorLineForce(self.turb, u_ph, self.simTime_id, self.dt, self.turb.index, self.mpi_u_fluid)

                # temp_uph = self.turb.build_actuator_lines(self.fs, mpi_u_fluid_ph, self.inflow_angle)
                # temp_umh = self.turb.build_actuator_lines(self.fs, mpi_u_fluid_mh, self.inflow_angle)

                self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_ph)
                temp_uph = self.turb.build_actuator_lines(self.fs, u_k, self.inflow_angle)

                self.turb.mpi_u_fluid_constant.assign(mpi_u_fluid_mh)
                temp_umh = self.turb.build_actuator_lines(self.fs, u_k, self.inflow_angle)

                self.turb.mpi_u_fluid_constant.assign(Constant(mpi_u_fluid_buff, name="temp_u_f"))         
                dtf_du = (temp_uph.vector().get_local()-temp_umh.vector().get_local())/(2.0*h_mag)

                prepared["u_local"].append(dtf_du)

        return prepared

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):

        ### Get the control type and turbine index ###
        name, turbine_number, segment_index = block_variable.tag
        # print("calculating: " + name + "_" + repr(segment_index))

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
