from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Create an Instance of the Options ###
options = windse.initialize("params3Dterrain.yaml")

### Generate Domain ###
dom = windse.ImportedDomain()
# dom = windse.BoxDomain()

### Generate Wind Farm ###
farm = windse.GridWindFarm(dom)
# farm = windse.RandomWindFarm(dom)

# farm.Plot(False)

### Warp the mesh and refine ###
# dom.Warp(200,0.75)
region = [[-1000,1500],[-1000,1000],[0,300]]
dom.Refine(1,region=region)
# dom.Save()


print(len(dom.mesh.coordinates()[:]))
print(len(farm.dom.mesh.coordinates()[:]))
print(print(farm.dom.mesh.hmin()))
# exit()

### Function Space ###
fs = windse.LinearFunctionSpace(dom)

print(fs.Q.dim())

### Setup Boundary Conditions ###
bc = windse.PowerInflow(dom,fs)

### Generate the problem ###
problem = windse.StabilizedProblem(dom,farm,fs,bc)

### Solve ###
solver = windse.SteadySolver(problem)
solver.Solve()

### Output Results ###
# solver.Save()

######
# control = windse.CreateAxialControl(farm.ma,farm)
# bounds = windse.CreateAxialBounds(farm.ma,farm)
control = windse.CreateLayoutControl(farm.mx,farm.my,farm)
bounds = windse.CreateLayoutBounds(farm.mx,farm.my,farm)

J=windse.PowerFunctional(problem.tf,solver.u_next)
# rf=ReducedFunctional(J,control)

# def iter_cb(m):
# 	# if MPI.rank(mpi_comm_world()) == 0:
# 	print("m = ")
# 	for mm in m:
# 		print("Constant("+ str(mm)+ "),")

# m_opt=minimize(rf, method="L-BFGS-B", options = {"disp": True}, bounds = bounds, callback = iter_cb)
# # m_opt=minimize(rf, method="SLSQP", options = {"disp": True}, bounds = bounds, callback = iter_cb)

# print([float(mm) for mm in m_opt])
# # farm.mx,farm.my=splitSolution(m_opt,farm.numturbs)
# farm.ma = m_opt
# solver.Solve()

mtest=[]
for i in range(farm.numturbs):
    mtest.append((farm.mx[i]))
    mtest.append((farm.my[i]))

h = [Constant(0.001)]*(2*farm.numturbs)  # the direction of the perturbation
Jhat = ReducedFunctional(J, control)  
conv_rate = taylor_test(Jhat, mtest, h)
print(conv_rate)

# ### Output Results ###
solver.Save()
