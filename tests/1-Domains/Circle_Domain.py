import windse
from dolfin import *

### Alias Parameters ###
params = windse.windse_parameters

### Set General Parameters ###
params["general"]["name"] = "Circle_Domain_Test"

### Set Box Parameters ###
radius = 1200
params["domain"]["type"]      = "circle"
params["domain"]["mesh_type"] = "mshr"
params["domain"]["radius"]    = radius
params["domain"]["center"]    = [0.0, 0.0]
params["domain"]["nt"]        = 100
params["domain"]["res"]       = 100

### Initialize Parameters using those set above ###
windse.initialize(None)

### Create the Domain Object ###
dom = windse.CircleDomain()

### Check if the object is as expected ###
dom.Save()

### Check the number of verts and cells ###
if dom.mesh.num_cells() != 3:
    raise ValueError("Box domain constructed with unexpected number of cells")
if dom.mesh.num_vertices() != 19865:
    raise ValueError("Box domain constructed with unexpected number of vertices")

### Calculate Volume ###
V = FiniteElement('Lagrange', dom.mesh.ufl_cell(), 1)
V = FunctionSpace(dom.mesh,V)
u = Function(V)
u.vector()[:] = 1.0
Volume = assemble(u*dx)
True_Volume = (pi*radius**2.0)

### Check the mesh volume ###
if abs((Volume-True_Volume)/True_Volume) > 1e-3:
    raise ValueError("Box domain constructed with unexpected volume")
