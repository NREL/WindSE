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
windse.initialize(params)

### Create the Domain Object ###
dom = windse.CircleDomain()

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

def test_volume():
    ### Calculate Volume ###
    Volume = assemble(u*dx)
    True_Volume = (pi*radius**2.0)

    ### Check the mesh volume ###
    if abs((Volume-True_Volume)/True_Volume) > 1e-3:
        print("Expected Volume: " + repr(True_Volume))
        print("Actual Volume:   " + repr(Volume))
        print("Percent Error:     " + repr(abs((Volume-True_Volume)/True_Volume)))
        raise ValueError("Box domain constructed with unexpected volume")

def test_boundary():
    ### Check the initial boundary integral
    Boundary_Value = CalculateInflowBoundary(dom,u)
    True_Boundary_Value = -2363.411257268333
    if abs((Boundary_Value-True_Boundary_Value)/True_Boundary_Value) > 1e-3:
        print("Expected inflow: " + repr(True_Boundary_Value))
        print("Actual inflow:   " + repr(Boundary_Value))
        raise ValueError("Initial inflow integral returned unexpected value (test value is hard coded)")

def test_rotated_boundary():
    ### Rotate the boundary by 1.337 rads
    dom.RecomputeBoundaryMarkers(1.337)

    ### Test Rotated Boundary integral ###
    Boundary_Value = CalculateInflowBoundary(dom,u)
    True_Boundary_Value = -2077.458137099408
    if abs((Boundary_Value-True_Boundary_Value)/True_Boundary_Value) > 1e-3:
        print("Expected inflow: " + repr(True_Boundary_Value))
        print("Actual inflow:   " + repr(Boundary_Value))
        # raise ValueError("Initial inflow integral returned unexpected value (test value is hard coded)")