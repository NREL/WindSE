from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Initialize WindSE ###
windse.initialize("params_small.yaml")

### Generate Domain ###
dom = windse.BoxDomain()

### Generate Wind Farm ###
# farm = windse.RandomWindFarm(dom)
farm = windse.GridWindFarm(dom)

# farm.Plot(False)

dom.Warp(250,0.75)
region = [farm.ex_x,farm.ex_y,[0,250]]
dom.Refine(1,region=region)
# dom.Save()

### Refine Around the Turbines
farm.RefineTurbines(num=1,radius_multiplyer=1.3) 

### Function Space ###
fs = windse.LinearFunctionSpace(dom)

### Setup Boundary Conditions ###
bc = windse.PowerInflow(dom,fs)

### Generate the problem ###
problem = windse.StabilizedProblem(dom,farm,fs,bc)

### Solve ###
solver = windse.SteadySolver(problem)
solver.Solve()

# ### Output Results ###
# solver.Save()
# solver.Plot()