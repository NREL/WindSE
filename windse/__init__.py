""" 
This is the init file for WindSE. It handle importing all the
submodules and initializing the parameters.
"""

from windse.ParameterManager import windse_parameters

def initialize(loc):
    """
    This function initializes all the submodules in WindSE.

    Args:
        loc (str): This string is the location of the .yaml parameters file.

    """

    windse_parameters.Load(loc)

    global BaseHeight, ReducedFunctional
    if windse_parameters["general"].get("dolfin_adjoint", False):
        from windse.dolfin_adjoint_helper import BaseHeight, ReducedFunctional
    else:
        from windse.helper_functions import BaseHeight
        
    global BoxDomain, CylinderDomain, RectangleDomain, ImportedDomain, InterpolatedCylinderDomain
    from windse.DomainManager import BoxDomain, CylinderDomain, RectangleDomain, ImportedDomain, InterpolatedCylinderDomain

    global GridWindFarm, RandomWindFarm, ImportedWindFarm
    from windse.WindFarmManager import GridWindFarm, RandomWindFarm, ImportedWindFarm

    global LinearFunctionSpace, TaylorHoodFunctionSpace2D
    from windse.FunctionSpaceManager import LinearFunctionSpace, TaylorHoodFunctionSpace2D

    global PowerInflow, UniformInflow, UniformInflow2D
    from windse.BoundaryManager import PowerInflow, UniformInflow, UniformInflow2D

    global StabilizedProblem, TaylorHoodProblem2D
    from windse.ProblemManager import StabilizedProblem, TaylorHoodProblem2D

    global SteadySolver, MultiAngleSolver
    from windse.SolverManager import SteadySolver, MultiAngleSolver

    global CreateAxialControl, CreateAxialBounds, CreateLayoutControl, CreateLayoutBounds, CreateYawControl, CreateYawBounds, SplitSolution, PowerFunctional
    from windse.OptimizationManager import CreateAxialControl, CreateAxialBounds, CreateLayoutControl, CreateLayoutBounds, CreateYawControl, CreateYawBounds, SplitSolution, PowerFunctional

