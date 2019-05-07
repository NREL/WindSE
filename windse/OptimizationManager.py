"""
The OptimizationManager submodule contains all the required function for
optimizing via dolfin-adjoint. To use dolfin-adjoin set::

    general: 
        dolfin_adjoint: True

in the param.yaml file.

Todo:
    * Read through an update the docstrings for these functions.
    * Create specific optimization classes.
"""

import __main__
import os

### Get the name of program importing this package ###
main_file = os.path.basename(__main__.__file__)

### This checks if we are just doing documentation ###
if main_file != "sphinx-build":
    from dolfin import *
    import numpy as np

    ### Import the cumulative parameters ###
    from windse import windse_parameters

    ### Check if we need dolfin_adjoint ###
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from dolfin_adjoint import *
    
def CreateLayoutControl(mx,my,farm):
    """
    Creates the controls from a list of values

    Args:
        m (list): a list of values to optimize.
    """
    m=[]
    for i in range(farm.numturbs):
        m.append(Control(mx[i]))
        m.append(Control(my[i]))

    return m

def CreateAxialControl(ma,farm):
    """
    Creates the controls from a list of values

    Args:
        m (list): a list of values to optimize.
    """
    m=[]
    for i in range(farm.numturbs):
        m.append(Control(ma[i]))
    return m
    
def CreateYawControl(myaw,farm):
    """
    Creates the controls from a list of values
    Args:
        m (list): a list of values to optimize.
    """
    m=[]
    for i in range(farm.numturbs):
        m.append(Control(myaw[i]))
    return m

def CreateAxialBounds(ma,farm):
    """
    Creates the optimization bounds for axial induction.

    Args:
        m (list): a list of controls
    """
    ub=[]
    lb=[]
    for i in range(farm.numturbs):
        lb.append(Constant(0))
        ub.append(Constant(1.))
        
    bounds = [lb,ub]
    return bounds

def CreateLayoutBounds(mx,my,farm):
        ub=[]
        lb=[]
        for i in range(farm.numturbs):
            lb.append(Constant((farm.ex_x[0] + farm.radius[i])))
            lb.append(Constant((farm.ex_y[0] + farm.radius[i])))
            ub.append(Constant((farm.ex_x[1] - farm.radius[i])))
            ub.append(Constant((farm.ex_y[1] - farm.radius[i])))
            
        bounds = [lb,ub]
        return bounds
        
def CreateYawBounds(myaw,farm):
    """
    Creates the optimization bounds for axial induction.

    Args:
        m (list): a list of controls
    """
    ub=[]
    lb=[]
    for i in range(farm.numturbs):
        lb.append(Constant(-pi/4.))
        ub.append(Constant(pi/4.))
        
    bounds = [lb,ub]
    return bounds

def SplitSolution(m_opt,numturbs):
    mx_opt = []
    my_opt = []
    j=0
    for i in range(numturbs):
        mx_opt.append(m_opt[j])
        j+=1
        my_opt.append(m_opt[j])
        j+=1
        
    return mx_opt,my_opt

def PowerFunctional(tf,u):
    """
    Creates the power functional that will be optimized

    Args:
        tf (dolfin.Function): Turbine Force function
        u (dolfin.Function): Velocity vector.
    """
    #how to handle rotation?
    # J=Functional(tf*u[0]**3*dx)
    J=assemble(-dot(tf,u)*dx)

    return J
