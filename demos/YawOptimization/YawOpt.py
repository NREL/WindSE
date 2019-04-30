from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Create an Instance of the Options ###
windse.initialize("params.yaml")

### Generate Domain ###
dom = windse.BoxDomain()

### Generate Wind Farm ###
farm = windse.ImportedWindFarm(dom)
# farm.Plot()

### Warp the mesh and refine ###
dom.Warp(200,0.75)
# dom.Save()

print(len(dom.mesh.coordinates()[:]))
print(len(farm.dom.mesh.coordinates()[:]))

### Function Space ###
fs = windse.LinearFunctionSpace(dom)

print(fs.W.dim())

### Setup Boundary Conditions ###
bc = windse.PowerInflow(dom,fs)

### Generate the problem ###
problem = windse.StabilizedProblem(dom,farm,fs,bc)

### Solve ###
solver = windse.SteadySolver(problem)
solver.Solve()

control = windse.CreateYawControl(farm.myaw, farm)
bounds = windse.CreateYawBounds(farm.ma, farm)
J=windse.PowerFunctional(problem.tf,solver.u_next)

rf=ReducedFunctional(J,control)
print(J)
print(float(J))
dJdma= compute_gradient(J, control, options={"newton_solver":{"linear_solver": "mumps"}})
print([float(dd) for dd in dJdma])

# def iter_cb(m):
# 	# if MPI.rank(mpi_comm_world()) == 0:
# 	print("m = ")
# 	for mm in m:
# 		print("Constant("+ str(mm)+ "),")

# m_opt=minimize(rf, method="L-BFGS-B", options = {"disp": True}, bounds = bounds, callback = iter_cb)
# print([float(mm) for mm in m_opt])
# farm.ma = m_opt
# solver.Solve()

# h = [Constant(0.001),Constant(0.001)]  # the direction of the perturbation
# Jhat = ReducedFunctional(J, control)  
# conv_rate = taylor_test(Jhat, farm.ma, h)
# print(conv_rate)

### Output Results ###
solver.Save()
print("finished")