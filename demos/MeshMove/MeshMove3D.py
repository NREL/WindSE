from dolfin import *
from dolfin_adjoint import *
import windse
import numpy as np
import copy

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt


parameters['form_compiler']['quadrature_degree'] = 6
set_log_level(15)

### Create an Instance of the Options ###
windse.initialize("params3D.yaml")

### Generate Simple Domain ###
dom = windse.BoxDomain()
dom.Save(val=0)

### Warp Mesh. This is always fine since the deltaZ doesn't change with x ###
dom.Warp(200,0.75)
dom.Save(val=1)

### Refine in the middle. This creates hanging nodes. ###
region = [[-100,100],[-100,100],[0,150]]
dom.Refine(1,region=region)
dom.Save(val=2)

# # dom.boundary_markers = MeshFunction("double", dom.mesh, dom.mesh.topology().dim() - 1)
# for facet in facets(dom.mesh):
#     ind = facet.index()
#     if dom.boundary_markers[ind] > 6:
#         dom.boundary_markers.set_value(ind,0)


### Save a copy of the mesh for returning to this state ###
x_temp = copy.deepcopy(dom.mesh.coordinates()[:,0])
y_temp = copy.deepcopy(dom.mesh.coordinates()[:,1])
z_temp = copy.deepcopy(dom.mesh.coordinates()[:,2])

### Define the transformation. It has very large gradients
def bed(x,y):
    return -100*np.sin(np.pi*x/250.0)*np.sin(np.pi*y/250.0)
    # return 150*np.sin(np.pi*x/400.0)
def transform(x,y,z,z0,z1):
    return (z1-bed(x,y))*(z-z0)/(z1-z0)+bed(x,y)

### Deform the mesh by directly changing the mesh coordinates
z_new = transform(x_temp,y_temp,z_temp,0.,500.)
dom.mesh.coordinates()[:,0]=x_temp
dom.mesh.coordinates()[:,1]=y_temp
dom.mesh.coordinates()[:,2]=z_new
dom.mesh.bounding_box_tree().build(dom.mesh)
dom.Save(val=3)

### Return to the original mesh ###
dom.mesh.coordinates()[:,0]=x_temp
dom.mesh.coordinates()[:,1]=y_temp
dom.mesh.coordinates()[:,2]=z_temp
dom.mesh.bounding_box_tree().build(dom.mesh)
dom.Save(val=4)

### Create a boundary mesh ###
b_mesh = BoundaryMesh(dom.mesh,"exterior")

### Move the boundary mesh ###
x_hd = copy.deepcopy(b_mesh.coordinates()[:,0])
y_hd = copy.deepcopy(b_mesh.coordinates()[:,1])
z_hd = copy.deepcopy(b_mesh.coordinates()[:,2])
z_hd = transform(x_hd,y_hd,z_hd,0.,500.)
b_mesh.coordinates()[:,0]=x_hd
b_mesh.coordinates()[:,1]=y_hd
b_mesh.coordinates()[:,2]=z_hd
b_mesh.bounding_box_tree().build(b_mesh)
plot(b_mesh)
plt.show()

### Move mesh using a hd boundary mesh ###
ALE.move(dom.mesh,b_mesh)
dom.Save(val=5)