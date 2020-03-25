import windse
from dolfin import *

### Alias Parameters ###
params = windse.windse_parameters

### Set General Parameters ###
params["general"]["name"] = "Box_Domain_Test"

### Set Box Parameters ###
Lx = 1000.0
Ly = 1000.0
Lz = 500
params["domain"]["type"]    = "box"
params["domain"]["x_range"] = [-Lx, Lx]
params["domain"]["y_range"] = [-Ly, Ly]
params["domain"]["z_range"] = [0.0, Lz]
params["domain"]["nx"]      = 24
params["domain"]["ny"]      = 24
params["domain"]["nz"]      = 6

### Initialize Parameters using those set above ###
windse.initialize(None)

### Create the Domain Object ###
dom = windse.BoxDomain()

### Check if the object is as expected ###
dom.Save()

### Check the number of verts and cells ###
if dom.mesh.num_cells() != 20736:
    raise ValueError("Box domain constructed with unexpected number of cells")
if dom.mesh.num_vertices() != 4375:
    raise ValueError("Box domain constructed with unexpected number of vertices")

### Calculate Volume ###
V = FiniteElement('Lagrange', dom.mesh.ufl_cell(), 1)
V = FunctionSpace(dom.mesh,V)
u = Function(V)
u.vector()[:] = 1.0
Volume = assemble(u*dx)
True_Volume = (2*Lx*2*Ly*Lz)

### Check the mesh volume ###
if abs((Volume-True_Volume)/True_Volume) > 1e-10:
    raise ValueError("Box domain constructed with unexpected volume")
