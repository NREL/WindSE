from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Create an Instance of the Options ###
windse.initialize("params.yaml")

### Generate Domain ###
dom = windse.RectangleDomain()

### Generate Wind Farm ###
farm = windse.GridWindFarm(dom)
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
solver.Save()
