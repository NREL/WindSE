import sys
import numpy as np
import time
import os.path as osp
import argparse
import sys
import windse

def run_driver(params_loc=None):
    # Run the model at least once
    params, problem, solver = run_model(params_loc=params_loc)

    ### Perform Optimization ###
    if params.get("optimization",{}):
        opt=windse.Optimizer(solver)
        if params["optimization"].get("taylor_test",False):
            opt.TaylorTest()

        if params["optimization"].get("optimize",True):
            opt.Optimize()

        if params["optimization"].get("gradient",True):
            opt.Gradient()
    else:
        raise NameError("paramters doesn't contain any optimization options")

def run_model(params_loc=None, comm=None):
    params = initialize_analysis(params_loc=params_loc, comm=comm)
    problem = setup_problem(params)
    solver = solve_problem(params, problem)

    return params, problem, solver

def initialize_analysis(params_loc=None, comm=None):

    # You need to parse the YAML file before initializing WindSE
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
    if comm is None:
        pass
    else:
        windse.windse_parameters['comm'] = comm

    params=windse.windse_parameters

    return params

def setup_problem(params):
    ### Setup the Domain ###
    if params["domain"].get("interpolated",False):
        dom_dict = {"box":windse.InterpolatedBoxDomain,
                    "cylinder":windse.InterpolatedCylinderDomain
        }
    else:
        dom_dict = {"box":windse.BoxDomain,
                    "rectangle":windse.RectangleDomain,
                    "cylinder":windse.CylinderDomain,
                    "circle":windse.CircleDomain,
                    "imported":windse.ImportedDomain}
    dom = dom_dict[params["domain"]["type"]]()

    ### Setup the Wind farm ###
    farm_dict = {"grid":windse.GridWindFarm,
                 "random":windse.RandomWindFarm,
                 "imported":windse.ImportedWindFarm}
    farm = farm_dict[params["wind_farm"]["type"]](dom)
    farm.Plot(params["wind_farm"].get("display",False))

    ### Move and refine the mesh
    if "refine" in params.keys():
        warp_type      = params["refine"].get("warp_type",None)
        warp_strength  = params["refine"].get("warp_strength",None)
        warp_height    = params["refine"].get("warp_height",None)
        warp_percent   = params["refine"].get("warp_percent",None)
        farm_num       = params["refine"].get("farm_num",0)
        farm_type      = params["refine"].get("farm_type","square")
        farm_factor    = params["refine"].get("farm_factor",1.0)
        farm_radius    = params["refine"].get("farm_radius",None)
        refine_custom  = params["refine"].get("refine_custom",None)
        turbine_num    = params["refine"].get("turbine_num",0)
        turbine_factor = params["refine"].get("turbine_factor",1.0)

        if warp_type == "smooth":
            dom.WarpSmooth(warp_strength)
        elif warp_type == "split":
            dom.WarpSplit(warp_height,warp_percent)


        if refine_custom is not None:
            for refine_data in refine_custom:
                if refine_data[1] == "full":
                    dom.Refine(refine_data[0])
                else:
                    region = farm.CalculateFarmRegion(refine_data[1],length=refine_data[2])
                    dom.Refine(refine_data[0],region=region,region_type=refine_data[1])

        if farm_num > 0:
            region = farm.CalculateFarmRegion(farm_type,farm_factor,length=farm_radius)
            dom.Refine(farm_num,region=region,region_type=farm_type)


        if turbine_num > 0:
            farm.RefineTurbines(turbine_num,turbine_factor)

    ### Finalize the Domain ###
    dom.Finalize()

    # dom.Save()
    # exit()
    ### Function Space ###
    func_dict = {"linear":windse.LinearFunctionSpace,
                 "taylor_hood":windse.TaylorHoodFunctionSpace}
    fs = func_dict[params["function_space"]["type"]](dom)

    # dom.Save()
    # exit()

    ### Setup Boundary Conditions ###
    bc_dict = {"uniform":windse.UniformInflow,
               "power":windse.PowerInflow,
               "log":windse.LogLayerInflow}
    bc = bc_dict[params["boundary_condition"]["vel_profile"]](dom,fs,farm)

    ### Generate the problem ###
    prob_dict = {"stabilized":windse.StabilizedProblem,
                 "taylor_hood":windse.TaylorHoodProblem}
    problem = prob_dict[params["problem"]["type"]](dom,farm,fs,bc)#,opt=opt)

    print('type', type(problem))
    return problem

def solve_problem(params, problem):
    solve_dict = {"steady":windse.SteadySolver,
                  "multiangle":windse.MultiAngleSolver,
                  "importedvelocity":windse.TimeSeriesSolver}
    solver = solve_dict[params["solver"]["type"]](problem)
    solver.Solve()
    return solver


if __name__ == '__main__':
    pass
