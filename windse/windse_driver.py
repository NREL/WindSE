import sys
import os.path as osp
import argparse

import numpy as np

import windse

ALL_ACTIONS = ("run", "cost")
help_msg = """
Available commands:

    run      run windse with a specified params file
    cost     run windse with specified eddy viscocity

Type windse <command> --help for usage help on a specific command.
For example, windse run --help will list all running options.
"""
### Print the help message ###
def print_usage():
    print("Usage: %s <command> <options> <arguments>"
          % osp.basename(sys.argv[0]))
    print(help_msg)

### Pop first argument, check it is a valid action. ###
def get_action():
    if len(sys.argv) == 1:
        action = 'run'
    elif len(sys.argv)==2:
        if sys.argv[1] in ALL_ACTIONS:
            action = sys.argv.pop(1)
        else:
            action = 'run'
    elif len(sys.argv)==3:
        if not sys.argv[1] in ALL_ACTIONS:
            print_usage()
            sys.exit(1)
        action = sys.argv.pop(1)
    else:
        print_usage()
        sys.exit(1)
    return action

### Run the driver ###
def run_action(params_loc=None):
    if params_loc is None:
        parser = argparse.ArgumentParser(usage="windse run [options] params", formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("params", nargs='?', help='path to yaml file containing the WindSE parameters')
        args = parser.parse_args()

        ### Check if parameters was provided ###
        if args.params is None:
            params_loc = "params.yaml"
            print("Using default parameter location: ./params.yaml")
        else:
            params_loc = args.params
            print("Using parameter location: "+params_loc)
    else:
        print("Using parameter location: "+params_loc)

    ### Initialize WindSE ###
    windse.initialize(params_loc)
    params=windse.windse_parameters

    ### Setup the Domain ###
    dom_dict = {"box":windse.BoxDomain,
                "rectangle":windse.RectangleDomain,
                "cylinder":windse.CylinderDomain,
                "imported":windse.ImportedDomain,
                "interpolated":windse.InterpolatedCylinderDomain}
    dom = dom_dict[params["domain"]["type"]]()

    ### Setup the Wind farm ###
    farm_dict = {"grid":windse.GridWindFarm,
                 "random":windse.RandomWindFarm,
                 "imported":windse.ImportedWindFarm}
    farm = farm_dict[params["wind_farm"]["type"]](dom)
    farm.Plot(False)

    ### Move and refine the mesh
    if "refine" in params.keys():
        warp_height    = params["refine"].get("warp_height",None)
        warp_percent   = params["refine"].get("warp_percent",None)
        farm_num       = params["refine"].get("farm_num",0)
        farm_type      = params["refine"].get("farm_type","square")
        farm_factor    = params["refine"].get("farm_factor",1.0)
        farm_custom    = params["refine"].get("farm_custom",None)
        turbine_num    = params["refine"].get("turbine_num",0)
        turbine_factor = params["refine"].get("turbine_factor",1.0)

        if warp_height is not None:
            dom.Warp(warp_height,warp_percent)

        if farm_num > 0:
            if farm_custom is not None:
                region = farm_custom
            else:
                region = farm.CalculateFarmRegion(farm_type,farm_factor)

            dom.Refine(farm_num,region=region,region_type=farm_type)

        if turbine_num > 0:
            farm.RefineTurbines(turbine_num,turbine_factor)

    ### Finalize the Domain ###
    dom.Finalize()

    ### Function Space ###
    func_dict = {"linear":windse.LinearFunctionSpace,
                 "taylor_hood":windse.TaylorHoodFunctionSpace}
    fs = func_dict[params["function_space"]["type"]](dom)

    ### Setup Boundary Conditions ###
    bc_dict = {"uniform":windse.UniformInflow,
               "power":windse.PowerInflow}
    bc = bc_dict[params["boundary_condition"]["vel_profile"]](dom,fs)

    ### Turbulence Model ###
    tm_dict = {"mixing_length":windse.MixingLengthTurbulenceModel,
               "specified_nuT":windse.SpecifiedNuT,
              }
    tm = tm_dict[params["physics"]["turbulence_model"]]()

    ### Generate the problem ###
    prob_dict = {"stabilized":windse.StabilizedProblem,
                 "taylor_hood":windse.TaylorHoodProblem}
    problem = prob_dict[params["problem"]["type"]](dom,farm,fs,bc,tm)

    ### Solve ###
    solve_dict = {"steady":windse.SteadySolver,
                  "multiangle":windse.MultiAngleSolver}
    solver = solve_dict[params["solver"]["type"]](problem)
    solver.Solve()


    # import dolfin_adjoint as da
    # tape = da.get_working_tape()
    # tape.visualise()
    # exit()

    ### Perform Optimization ###
    if "optimization" in params.keys():
        opt=windse.Optimizer(solver)
        if params["optimization"].get("taylor_test",False):
            opt.TaylorTest()

        if params["optimization"].get("optimize",True):
            opt.Optimize()


def main():
    actions = {"run": run_action,}
    actions[get_action()]()

if __name__ == "__main__":
    main()
