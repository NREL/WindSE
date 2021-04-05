import numpy as np
import time
import os.path as osp
import argparse
import sys
import windse

def DefaultParameters():
    """
    return the default parameters list
    """
    return windse.windse_parameters.defaults

def BlankParameters():
    """
    returns a nested dictionary that matches the first level of the parameters dictionary
    """
    params = {}
    params["general"] = {}
    params["domain"] = {}
    params["wind_farm"] = {}
    params["function_space"] = {}
    params["boundary_conditions"] = {}
    params["problem"] = {}
    params["solver"] = {}
    params["optimization"] = {}
    return params


def Initialize(params_loc=None):
    """
    This function initialized the windse parameters.

    Parameters
    ----------
        params_loc : str
            the location of the parameter yaml file.

    Returns
    -------
        params : windse.Parameters
            an overloaded dict containing all parameters.
    """
    parser = argparse.ArgumentParser(usage="windse run [options] params", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("params", nargs='?', help='path to yaml file containing the WindSE parameters')
    parser.add_argument('-p', dest='updated_parameters', action='append',default=[], help='use this to override a parameter in the yaml file')
    args, unknown = parser.parse_known_args()

    if params_loc is None:
        ### Check if parameters was provided ###
        if args.params is None:
            params_loc = "params.yaml"
            print("Using default parameter location: ./params.yaml")
        else:
            params_loc = args.params
            print("Using parameter location: "+params_loc)

    ### Initialize WindSE ###
    windse.initialize(params_loc,updated_parameters=args.updated_parameters)
    
    params=windse.windse_parameters

    return params

def BuildDomain(params):
    """
    This function build the domain and wind farm objects.

    Parameters
    ----------
    params : ``windse.Parameters``
        an overloaded dict containing all parameters.

    Returns
    -------
    dom : ``windse.GenericDomain``
        the domain object that contains all mesh related information.
    farm : ``windse.GenericWindFarm`` 
        the wind farm object that contains the turbine information.
    """

    ### Build Domain ###
    if params["domain"]["interpolated"]:
        dom_dict = {"box":windse.InterpolatedBoxDomain,
                    "cylinder":windse.InterpolatedCylinderDomain,
                    "imported":windse.ImportedDomain
        }
    else:
        dom_dict = {"box":windse.BoxDomain,
                    "rectangle":windse.RectangleDomain,
                    "cylinder":windse.CylinderDomain,
                    "circle":windse.CircleDomain,
                    "imported":windse.ImportedDomain}
    dom = dom_dict[params["domain"]["type"]]()


    #### Build Farm
    farm_dict = {"grid":windse.GridWindFarm,
                 "random":windse.RandomWindFarm,
                 "imported":windse.ImportedWindFarm,
                 "empty":windse.EmptyWindFarm}
    farm = farm_dict[params["wind_farm"]["type"]](dom)

    if dom.type != "imported":
        ### warp and refine the mesh
        windse.WarpMesh(dom)
        windse.RefineMesh(dom,farm)

        ### Finalize the Domain ###
        dom.Finalize()

    return dom, farm

def BuildProblem(params,dom,farm):
    """
    This function compiles everything into a single problem object and build the variational problem functional.

    Args:
        params (windse.Parameters): an overloaded dict containing all parameters.
        dom (windse.GenericDomain): the domain object that contains all mesh related information.
        farm (windse.GenericWindFarm): the wind farm object that contains the turbine information.

    Returns:
        problem (windse.GenericProblem): contains all information about the simulation.
    """
    func_dict = {"linear":windse.LinearFunctionSpace,
                 "taylor_hood":windse.TaylorHoodFunctionSpace}
    fs = func_dict[params["function_space"]["type"]](dom)

    # dom.Save()
    # exit()
    
    ### Setup Boundary Conditions ###
    bc_dict = {"uniform":windse.UniformInflow,
               "power":windse.PowerInflow,
               "log":windse.LogLayerInflow,
               "turbsim":windse.TurbSimInflow}
    bc = bc_dict[params["boundary_conditions"]["vel_profile"]](dom,fs,farm)

    ### Generate the problem ###
    prob_dict = {"stabilized":windse.StabilizedProblem,
                 "steady":windse.StabilizedProblem,
                 "taylor_hood":windse.TaylorHoodProblem,
                 "iterative_steady":windse.IterativeSteady,
                 "unsteady":windse.UnsteadyProblem}
    problem = prob_dict[params["problem"]["type"]](dom,farm,fs,bc)#,opt=opt)

    return problem

def BuildSolver(params,problem):
    """
    This function builds the solver object. Solve with solver.Solve()

    Parameters
    ----------
        params : windse.Parameters
            an overloaded dict containing all parameters.
        problem : windse.GenericProblem
            contains all information about the simulation.

    Returns
    -------
        solver : windse.GenericSolver
            contains the solver routines.
    """

    if isinstance(params["boundary_conditions"]["inflow_angle"],list):
        params["solver"]["type"] = "multiangle"
    solve_dict = {"steady":windse.SteadySolver,
                  "iterative_steady":windse.IterativeSteadySolver,
                  "unsteady":windse.UnsteadySolver,
                  "multiangle":windse.MultiAngleSolver,
                  "imported_inflow":windse.TimeSeriesSolver}
    solver = solve_dict[params["solver"]["type"]](problem)

    return solver

def SetupSimulation(params_loc=None):
    """
    This function automatically sets up the entire simulation. Solve with solver.Solve()

    Parameters
    ----------
        params_loc : str
            the location of the parameter yaml file.

    Returns
    -------
        params : windse.Parameters
            an overloaded dict containing all parameters.
        problem : windse.GenericProblem
            contains all information about the simulation.
        solver : windse.GenericSolver
            contains the solver routines. Solve with solver.Solve()
    """
    params = Initialize(params_loc)
    dom, farm = BuildDomain(params)
    problem = BuildProblem(params,dom,farm)
    solver = BuildSolver(params,problem)

    farm.PlotFarm(params["wind_farm"]["display"])
    if farm.chord is not None:
        farm.PlotChord(params["wind_farm"]["display"])


    return params, problem, solver
