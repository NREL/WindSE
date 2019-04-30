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
windse.initialize("params.yaml")

### Generate Simple Domain ###
dom = windse.RectangleDomain()
dom.Save(val=0)

### Warp Mesh. This is always fine since the deltaZ doesn't change with x ###
dom.Warp(200,0.75)
dom.Save(val=1)

### Refine in the middle. This creates hanging nodes. ###
region = [[-100,100],[0,150],[0,150]]
dom.Refine(1,region=region)
dom.Save(val=2)


# print(dir(dom.boundary_markers))
# print(type(dom.boundary_markers.where_equal(2)))
# print(dom.boundary_markers.where_equal(2))
# print(dir(Facet(dom.mesh,33)))
# ids = Facet(dom.mesh,33).entities(0)
# coord = dom.mesh.coordinates()[ids]
# print(coord[0,0],coord[1,1])
# if all(coord[:,1]<10):
#     print('hooray')
# exit()

### Recolor boundary function to be the original z values ###
dom.boundary_markers = MeshFunction("double", dom.mesh, dom.mesh.topology().dim() - 1)
for facet in facets(dom.mesh):
    ind = facet.index()
    dom.boundary_markers.set_value(ind,facet.midpoint()[1])
dom.Save(val=3)

### Save a copy of the mesh for returning to this state ###
x_temp = copy.deepcopy(dom.mesh.coordinates()[:,0])
z_temp = copy.deepcopy(dom.mesh.coordinates()[:,1])

### Define the transformation. It has very large gradients
def bed(x):
    return -100*np.sin(np.pi*x/200.0)
    # return 150*np.sin(np.pi*x/400.0)
def transform(x,z,z0,z1):
    return (z1-bed(x))*(z-z0)/(z1-z0)+bed(x)

### Deform the mesh by directly changing the mesh coordinates
z_new = transform(x_temp,z_temp,0.,500.)
dom.mesh.coordinates()[:,0]=x_temp
dom.mesh.coordinates()[:,1]=z_new
dom.mesh.bounding_box_tree().build(dom.mesh)
dom.Save(val=4)

### Return to the original mesh ###
dom.mesh.coordinates()[:,0]=x_temp
dom.mesh.coordinates()[:,1]=z_temp
dom.mesh.bounding_box_tree().build(dom.mesh)
dom.Save(val=5)

### This will create a displacement vector for every node ###
fs = windse.TaylorHoodFunctionSpace2D(dom)
displacement_x = Function(fs.V0)
displacement_z = Function(fs.V1)
dof_coords = fs.V0.tabulate_dof_coordinates()
dis_temp_x = np.zeros(len(displacement_x.vector()[:]))
dis_temp_z = np.zeros(len(displacement_z.vector()[:]))
for i in range(len(dof_coords)):
    dis_temp_x[i] = 0.0
    dis_temp_z[i] = transform(dof_coords[i,0],dof_coords[i,1],0.,500.)-dof_coords[i,1]
displacement_x.vector()[:]=dis_temp_x
displacement_z.vector()[:]=dis_temp_z
displacement=Function(fs.V)
fs.VelocityAssigner.assign(displacement,[displacement_x,displacement_z])

### Plot this displacement vector field ###
plot(displacement)
plt.show()

### Move the mesh using ALE ###
ALE.move(dom.mesh,displacement)
dom.Save(val=6)

### Return to original mesh ###
dom.mesh.coordinates()[:,0]=x_temp
dom.mesh.coordinates()[:,1]=z_temp
dom.mesh.bounding_box_tree().build(dom.mesh)
dom.Save(val=7)

### Create a boundary mesh ###
b_mesh = BoundaryMesh(dom.mesh,"exterior")

### Move the boundary mesh ###
x_hd = copy.deepcopy(b_mesh.coordinates()[:,0])
z_hd = copy.deepcopy(b_mesh.coordinates()[:,1])
z_hd = transform(x_hd,z_hd,0.,500.)
b_mesh.coordinates()[:,0]=x_hd
b_mesh.coordinates()[:,1]=z_hd
b_mesh.bounding_box_tree().build(b_mesh)
plot(b_mesh)
plt.show()

### Move mesh using a hd boundary mesh ###
ALE.move(dom.mesh,b_mesh)
dom.Save(val=8)