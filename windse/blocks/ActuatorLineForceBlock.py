from pyadjoint.block import Block

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
            adj_output = np.float64(np.inner(comp_vec, adj_vec))

            recv_buff = np.zeros(1, dtype=np.float64)
            self.problem.params.comm.Allreduce(adj_output, recv_buff)
            adj_output = recv_buff

            # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
            return np.array(adj_output)
        else:
            comp_vec = prepared[name][:, segment_index]
            adj_vec = adj_inputs[0].get_local()
            adj_output = np.float64(np.inner(comp_vec, adj_vec))

            recv_buff = np.zeros(1, dtype=np.float64)
            self.problem.params.comm.Allreduce(adj_output, recv_buff)
            adj_output = recv_buff

            # adj_output = dolfin.MPI.sum(dolfin.MPI.comm_world,adj_output)
            return np.array(adj_output)