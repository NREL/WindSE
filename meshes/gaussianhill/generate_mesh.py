from dolfin import *
import numpy as np
from scipy.interpolate import interp2d,interp1d, PchipInterpolator

### This import improves the plotter functionality on Mac ###
import matplotlib
matplotlib.use('TKAgg')

### Import matplotlib for plots ###
import matplotlib.pyplot as plt
parameters["refinement_algorithm"] = "plaza_with_parent_facets"

### Generate the meshes that will be manipulated
Lx0 = -2500
Lx1 =  2500
Ly0 = -2500
Ly1 =  2500
Lz0 =  0
Lz1 =  800
mesh = BoxMesh(Point(Lx0,Ly0,Lz0),Point(Lx1,Ly1,Lz1),30,30,10)

### Create Sample Function space to calculated dofs ### 
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, P2*P1)
print(W.dim())

### Load mesh grid data for the hill ###
print("Load Nodes")
hill_data = np.loadtxt("GaussianMeshNodes.txt")
x_data = hill_data[:,0]
y_data = hill_data[:,1]
z_data = hill_data[:,2]

### Generate the interpolation function for the hill ###
print("Build Interpolation Function")
mesh_interp = interp2d(x_data,y_data,z_data,'cubic')
print(mesh_interp(0,0))
print("Finished Interpolation Function")

### Generate a cubic spline for the density function ###
# cubic_spline = PchipInterpolator([0,0.7,0.8,1],[0,0.25,0.35,1])
cubic_spline = interp1d([0,0.7,0.8,1],[0,0.25,0.35,1])
# cubic_spline = PchipInterpolator([0,1],[0,1])

### Plot the spline to make sure it looks good ##
# x = np.linspace(0,1,50)
# y = cubic_spline(x)
# plt.plot(x,y)
# plt.show()

### Define an spline function distribution ###
def spline_dist_func(z,top,bottom,power):
    return (top-bottom)*cubic_spline((z-Lz0)/(Lz1-Lz0))**power+bottom

### Define a method for moving the mesh ###
def mesh_convert(x,y,z,method,power):
    n = len(x)
    x_new = np.zeros(n)
    y_new = np.zeros(n)
    z_new = np.zeros(n)
    for i in range(n):
        x_new[i] = x[i]
        y_new[i] = y[i]
        z_new[i] = method(z[i],Lz1,mesh_interp(x[i],y[i]),power)
    return [x_new,y_new,z_new]

### Move Mesh based on a Spline distribution ###
print("Move Mesh (Spline)")
x_mesh, y_mesh, z_mesh = mesh.coordinates()[:, 0], mesh.coordinates()[:, 1], mesh.coordinates()[:, 2]
x_bar, y_bar, z_bar = mesh_convert(x_mesh,y_mesh,z_mesh,spline_dist_func,1.0)
new_coords = np.array([x_bar,y_bar,z_bar]).transpose()
mesh.coordinates()[:] = new_coords

### Save mesh ###
print("Save Mesh (Spline)")

### Define the boundary subdomains ###
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2],Lz1) and on_boundary
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], mesh_interp(x[0],x[1])[0]) and on_boundary
class Front(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],Lx0) and on_boundary
class Back(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0],Lx1) and on_boundary
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1],Ly0) and on_boundary
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1],Ly1) and on_boundary
topb = Top()
bottomb = Bottom()
frontb = Front()
backb = Back()
leftb = Left()
rightb = Right()

### Create facet function for boundaries ###
boundaries = MeshFunction("size_t",mesh,mesh.geometry().dim()-1)
boundaries.set_all(0)
topb.mark(boundaries,1,check_midpoint=False)
bottomb.mark(boundaries,2,check_midpoint=False)
frontb.mark(boundaries,3,check_midpoint=False)
backb.mark(boundaries,4,check_midpoint=False)
leftb.mark(boundaries,5,check_midpoint=False)
rightb.mark(boundaries,6,check_midpoint=False)

### Save Mesh andBoundary facet function ###
vfile = File("pvd/mesh.pvd")
vfile << mesh

vfile = File("pvd/boundaries.pvd")
vfile << boundaries

vfile = File("xml/boundaries.xml.gz")
vfile << boundaries

vfile = File("xml/mesh.xml.gz")
vfile << mesh

