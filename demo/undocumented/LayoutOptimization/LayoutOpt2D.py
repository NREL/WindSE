from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Create an Instance of the Options ###
windse.initialize("params2Dsmall.yaml")

### Generate Domain ###
dom = windse.RectangleDomain()

### Generate Wind Farm ###
farm = windse.RandomWindFarm(dom)
# farm.Plot()

### Warp the mesh and refine ###
# dom.Warp(80,80,0.75)
# region = [[farm.ex_x[0],dom.x_range[1]],farm.ex_y,farm.ex_z]
# dom.Refine(1,local=True)
# dom.Save()

### Function Space ###
fs = windse.TaylorHoodFunctionSpace2D(dom)

### Setup Boundary Conditions ###
bc = windse.UniformInflow2D(dom,fs)

### Generate the problem ###
problem = windse.TaylorHoodProblem2D(dom,farm,fs,bc)

### Solve ###
solver = windse.SteadySolver(problem)
solver.Solve()

### Output Results ###
# solver.Save()
# exit()

####

# control = windse.CreateAxialControl(farm.ma,farm)
# bounds = windse.CreateAxialBounds(farm.ma,farm)
control = windse.CreateLayoutControl(farm.mx,farm.my,farm)
bounds = windse.CreateLayoutBounds(farm.mx,farm.my,farm)

J=windse.PowerFunctional(problem.tf,solver.u_next)
rf=ReducedFunctional(J,control)
# # # print(J)
# # # print(float(J))
# # # dJdma= compute_gradient(J, control)
# # # print([float(dd) for dd in dJdma])

def iter_cb(m):
	# if MPI.rank(mpi_comm_world()) == 0:
	print("m = ")
	for mm in m:
		print("Constant("+ str(mm)+ "),")

m_opt=minimize(rf, method="L-BFGS-B", options = {"disp": True}, bounds = bounds, callback = iter_cb)
# m_opt=minimize(rf, method="SLSQP", options = {"disp": True}, bounds = bounds, callback = iter_cb)

print([float(mm) for mm in m_opt])
farm.mx,farm.my=splitSolution(m_opt,farm.numturbs)
# farm.ma = m_opt
solver.Solve()

# h = [Constant(0.001)]*farm.numturbs  # the direction of the perturbation
# Jhat = ReducedFunctional(J, control)  
# conv_rate = taylor_test(Jhat, farm.ma, h)
# print(conv_rate)

# ### Output Results ###
solver.Save()
