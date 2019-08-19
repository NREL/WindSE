""" 
This is the init file for WindSE. It handle importing all the
submodules and initializing the parameters.
"""

from windse.ParameterManager import windse_parameters
def initialize(loc,updated_parameters=[]):
    """
    This function initializes all the submodules in WindSE.

    Args:
        loc (str): This string is the location of the .yaml parameters file.

    """

    windse_parameters.Load(loc,updated_parameters=updated_parameters)

    global BaseHeight, Optimizer#, ReducedFunctional
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from windse.dolfin_adjoint_helper import BaseHeight#, ReducedFunctional
        from windse.OptimizationManager import Optimizer
    else:
        from windse.helper_functions import BaseHeight
        
    global BoxDomain, CylinderDomain, CircleDomain, RectangleDomain, ImportedDomain, InterpolatedCylinderDomain, InterpolatedBoxDomain
    from windse.DomainManager import BoxDomain, CylinderDomain, CircleDomain, RectangleDomain, ImportedDomain, InterpolatedCylinderDomain, InterpolatedBoxDomain

    global GridWindFarm, RandomWindFarm, ImportedWindFarm
    from windse.WindFarmManager import GridWindFarm, RandomWindFarm, ImportedWindFarm

    global LinearFunctionSpace, TaylorHoodFunctionSpace
    from windse.FunctionSpaceManager import LinearFunctionSpace, TaylorHoodFunctionSpace

    global PowerInflow, UniformInflow, LogLayerInflow
    from windse.BoundaryManager import PowerInflow, UniformInflow, LogLayerInflow

    global StabilizedProblem, TaylorHoodProblem
    from windse.ProblemManager import StabilizedProblem, TaylorHoodProblem

    global SteadySolver, MultiAngleSolver
    from windse.SolverManager import SteadySolver, MultiAngleSolver

