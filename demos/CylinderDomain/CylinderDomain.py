from dolfin import *
import windse
import numpy as np

parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(20)

### Create an Instance of the Options ###
options = windse.initialize("params.yaml")

### Generate Domain ###
dom = windse.CylinderDomain()
# dom.Save()
# exit()



# exit()

### Generate Wind Farm ###
farm = windse.GridWindFarm(dom)
farm.Plot(False)

### Warp the mesh and refine ###
dom.Warp(200,0.75)
print(np.unique(dom.mesh.coordinates()[:,2]))
# region = [[-1500,1500],[-1500,1500],[0,300]]
region = [[1500],[0.0,0.0],[0,300]]
dom.Refine(1,region=region,region_type="circle")

### Refine Around the Turbines
farm.RefineTurbines(num_refinements=1,radius_multiplyer=1.4) 

### Function Space ###
fs = windse.LinearFunctionSpace(dom)

# print(fs.W.dim())
# dom.Save()
# exit()

### Setup Boundary Conditions ###
bc = windse.PowerInflow(dom,fs)

### Generate the problem ###
problem = windse.StabilizedProblem(dom,farm,fs,bc)

### Solve ###
# angles = np.linspace(0,np.pi/2,2)
# solver = windse.MultiAngleSolver(problem,angles)
# solver.Solve()

solver = windse.SteadySolver(problem)
solver.Solve()

# # exit()
# solver.Solve(iter_val=0.0)

# for i, theta in enumerate(angles):
#     if i > 0:
#         print("Computing for Wind Angle: "+repr(theta))
#         print("Adjusting for Angle "+repr(i+1)+" of "+repr(len(angles)))
#         solver.ChangeWindAngle(theta)
#         solver.Solve(iter_val=theta)
#         print("Finished Angle "+repr(i+1)+" of "+repr(len(angles)))

