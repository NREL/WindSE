from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np
import inspect

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Create an Instance of the Options ###
options = windse.initialize("GOparams.yaml")

### Generate Domain ###
dom = windse.ImportedDomain()

### Generate Wind Farm ###
farm = windse.GridWindFarm(dom)
farm.Plot(False)

### Warp the mesh and refine ###
region = [[-1000,1500],[-1000,1000],[0,300]]
# dom.Refine(1)
dom.Refine(1,region=region)
dom.Save()

### Function Space ###
fs = windse.LinearFunctionSpace(dom)

### Setup Boundary Conditions ###
bc = windse.PowerInflow(dom,fs)

### Generate the problem ###
problem = windse.StabilizedProblem(dom,farm,fs,bc)

### Solve ###
solver = windse.SteadySolver(problem)
solver.Solve()
solver.Save(val=0)

###### Optimization ######
control = windse.CreateLayoutControl(farm.mx,farm.my,farm)
bounds = windse.CreateLayoutBounds(farm.mx,farm.my,farm)

J=windse.PowerFunctional(problem.tf,solver.u_next)
rf=windse.ReducedFunctional(J,control)

def iter_cb(m):
    rf.iter_complete = True


m_opt=minimize(rf, method="L-BFGS-B", options = {"disp": True}, bounds = bounds, callback = iter_cb)
# m_opt=minimize(rf, method="SLSQP", options = {"disp": True}, bounds = bounds, callback = iter_cb)

print([float(mm) for mm in m_opt])
farm.mx,farm.my=splitSolution(m_opt,farm.numturbs)
# farm.ma = m_opt
solver.Solve()
solver.Save(val=solver.step)



###### Taylor Test ######
# control = windse.CreateLayoutControl(farm.mx,farm.my,farm)
# bounds = windse.CreateLayoutBounds(farm.mx,farm.my,farm)

# J=windse.PowerFunctional(problem.tf,solver.u_next)

# mtest=[]
# for i in range(farm.numturbs):
#     mtest.append((farm.mx[i]))
#     mtest.append((farm.my[i]))

# h = [Constant(0.001)]*(2*farm.numturbs)  # the direction of the perturbation
# Jhat = ReducedFunctional(J, control)  
# conv_rate = taylor_test(Jhat, mtest, h)
# print(conv_rate)










# ### SIMPLE Taylor testing code ###
# mesh = UnitSquareMesh(10, 10)
# V = FunctionSpace(mesh, 'CG', 1)

# def ground(x,y):
#     return float(x)*float(y)

# x = Constant(0.1,name="x")
# y = Constant(0.2,name="y")
# z = windse.BaseHeight(x,y,ground)

# unit = project(Constant(1.0),V)

# J = assemble((x+y+z)*unit*dx)

# h = [Constant(0.001),Constant(0.001)] 

# taylor_test(ReducedFunctional(J, [Control(x),Control(y)]), [x,y], h)