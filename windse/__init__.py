""" 
This is the init file for WindSE. It handle importing all the
submodules and initializing the parameters.
"""

import os
import __main__

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

from windse.ParameterManager import windse_parameters

def initialize(loc,updated_parameters=[]):
    """
    This function initializes all the submodules in WindSE.

    Args:
        loc (str): This string is the location of the .yaml parameters file.

    """

    windse_parameters.Load(loc,updated_parameters=updated_parameters)

    global BaseHeight, CalculateDiskTurbineForces, UpdateActuatorLineForce, RadialChordForce, Optimizer#, ReducedFunctional
    if windse_parameters["general"].get("dolfin_adjoint", False) or main_file in ["sphinx-build", "__main__.py"]:
        from windse.dolfin_adjoint_helper import BaseHeight, CalculateDiskTurbineForces, UpdateActuatorLineForce, RadialChordForce#, ReducedFunctional
        from windse.OptimizationManager import Optimizer
    else:
        from windse.helper_functions import BaseHeight, CalculateDiskTurbineForces, UpdateActuatorLineForce, RadialChordForce
    
    global BoxDomain, CylinderDomain, CircleDomain, RectangleDomain, ImportedDomain, InterpolatedCylinderDomain, InterpolatedBoxDomain, PeriodicDomain
    from windse.DomainManager import BoxDomain, CylinderDomain, CircleDomain, RectangleDomain, ImportedDomain, InterpolatedCylinderDomain, InterpolatedBoxDomain, PeriodicDomain

    global GridWindFarm, RandomWindFarm, ImportedWindFarm, EmptyWindFarm
    from windse.WindFarmManager import GridWindFarm, RandomWindFarm, ImportedWindFarm, EmptyWindFarm

    global RefineMesh, WarpMesh
    from windse.RefinementManager import RefineMesh, WarpMesh

    global LinearFunctionSpace, TaylorHoodFunctionSpace
    from windse.FunctionSpaceManager import LinearFunctionSpace, TaylorHoodFunctionSpace

    global PowerInflow, UniformInflow, LogLayerInflow, TurbSimInflow
    from windse.BoundaryManager import PowerInflow, UniformInflow, LogLayerInflow, TurbSimInflow

    global StabilizedProblem, TaylorHoodProblem, IterativeSteady, UnsteadyProblem
    from windse.ProblemManager import StabilizedProblem, TaylorHoodProblem, IterativeSteady, UnsteadyProblem

    global SteadySolver, IterativeSteadySolver, UnsteadySolver, MultiAngleSolver, TimeSeriesSolver
    from windse.SolverManager import SteadySolver, IterativeSteadySolver, UnsteadySolver, MultiAngleSolver, TimeSeriesSolver

