import windse
import pytest
import numpy as np
from dolfin import *
import windse_driver.driver_functions as df

###############################################################
######################## Setup Objects ########################
###############################################################

### Alias Parameters ###
params = df.BlankParameters()

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
windse.initialize(params)

### Create the Domain Object ###
dom = windse.BoxDomain()

### Check if the object is as expected ###
dom.Save()

### Create unit for integration ###
V = FiniteElement('Lagrange', dom.mesh.ufl_cell(), 1)
V = FunctionSpace(dom.mesh,V)
u = Function(V)
u.vector()[:] = 1.0

### Calculate inflow integral ###
def CalculateInflowBoundary(dom,u):
    xy = Expression("(x[0]+pi)*(x[1]+pi)",degree=2,pi=pi)
    ds = Measure('ds', subdomain_data=dom.boundary_markers)
    val = 0
    unique_id = np.unique(dom.boundary_markers.array())
    for blabel in dom.boundary_types["inflow"]: 
        if dom.boundary_names[blabel] in unique_id:
            val += assemble(xy*u*ds(dom.boundary_names[blabel]))/assemble(u*ds(dom.boundary_names[blabel]))
    return val







###############################################################
######################## Define Tests #########################
###############################################################

def test_mesh_values():
    True_Cells = 20736
    True_Verts = 4375
    if dom.mesh.num_cells() != True_Cells:
        print("Expected Cells: " + repr(True_Cells))
        print("Actual Cells:   " + repr(dom.mesh.num_cells()))
        raise ValueError("Box domain constructed with unexpected number of cells")
    if dom.mesh.num_vertices() != True_Verts:
        print("Expected x-values: " + repr(True_Verts))
        print("Actual x-values:   " + repr(dom.mesh.num_vertices()))
        raise ValueError("Box domain constructed with unexpected number of vertices")

def test_volume():
    ### Calculate Volume ###
    Volume = assemble(u*dx)
    True_Volume = (2*Lx*2*Ly*Lz)

    ### Check the mesh volume ###
    if abs((Volume-True_Volume)/True_Volume) > 1e-5:
        print("Expected Volume: " + repr(True_Volume))
        print("Actual Volume:   " + repr(Volume))
        print("Percent Error:     " + repr(abs((Volume-True_Volume)/True_Volume)))
        raise ValueError("Box domain constructed with unexpected volume")

def test_boundary():
    ### Check the initial boundary integral
    Boundary_Value = CalculateInflowBoundary(dom,u)
    True_Boundary_Value = -3111.983840386548
    if abs((Boundary_Value-True_Boundary_Value)/True_Boundary_Value) > 1e-5:
        print("Expected inflow: " + repr(True_Boundary_Value))
        print("Actual inflow:   " + repr(Boundary_Value))
        print("Percent Error:   " + repr(abs((Boundary_Value-True_Boundary_Value)/True_Boundary_Value)))
        raise ValueError("Initial inflow integral returned unexpected value (test value is hard coded)")

def test_rotated_boundary():
    ### Rotate the boundary by 1.337 rads
    dom.RecomputeBoundaryMarkers(1.337)

    ### Test Rotated Boundary integral ###
    Boundary_Value = CalculateInflowBoundary(dom,u)
    True_Boundary_Value = -3131.723049188729
    if abs((Boundary_Value-True_Boundary_Value)/True_Boundary_Value) > 1e-5:
        print("Expected inflow: " + repr(True_Boundary_Value))
        print("Actual inflow:   " + repr(Boundary_Value))
        print("Percent Error:   " + repr(abs((Boundary_Value-True_Boundary_Value)/True_Boundary_Value)))
        raise ValueError("Initial inflow integral returned unexpected value (test value is hard coded)")
