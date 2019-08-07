import sys
import numpy as np
import time
import os.path as osp
import argparse

ALL_ACTIONS = ("run")
help_msg = """
Available commands:

    run      run windse with a specified params file

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
    if len(sys.argv) <= 1:
        print_usage()
        sys.exit(1)
    if not sys.argv[1] in ALL_ACTIONS:
        if sys.argv[1] == "--help":
            print_usage()
            exit()
        else:
            print_usage()
            sys.exit(1)
    return sys.argv.pop(1)

### Run the driver ###
def run_action():
    tick = time.time()
    import windse
    
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

    ### Initialize WindSE ###
    windse.initialize(params_loc)
    
    params=windse.windse_parameters

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
        warp_height    = params["refine"].get("warp_height",None)
        warp_percent   = params["refine"].get("warp_percent",None)
        farm_num       = params["refine"].get("farm_num",0)
        farm_type      = params["refine"].get("farm_type","square")
        farm_factor    = params["refine"].get("farm_factor",1.0)
        farm_radius    = params["refine"].get("farm_radius",None)
        farm_custom    = params["refine"].get("farm_custom",None)
        turbine_num    = params["refine"].get("turbine_num",0)
        turbine_factor = params["refine"].get("turbine_factor",1.0)


        if farm_custom is not None:
            for refine_data in farm_custom:
                if refine_data[1] == "full":
                    dom.Refine(refine_data[0])
                else:
                    region = farm.CalculateFarmRegion(refine_data[1],length=refine_data[2])
                    dom.Refine(refine_data[0],region=region,region_type=refine_data[1])
        else:
            region = farm.CalculateFarmRegion(farm_type,farm_factor,length=farm_radius)
            dom.Refine(farm_num,region=region,region_type=farm_type)

###############################################
###############################################
###############################################
###############################################
        if warp_height is not None:
            dom.WarpNonlinear(1.2)
            # dom.Warp(warp_height,warp_percent)
###############################################
###############################################
###############################################
###############################################

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
               "power":windse.PowerInflow,
               "log":windse.LogLayerInflow}
    bc = bc_dict[params["boundary_condition"]["vel_profile"]](dom,fs)

    ### Generate the problem ###
    prob_dict = {"stabilized":windse.StabilizedProblem,
                 "taylor_hood":windse.TaylorHoodProblem}
    problem = prob_dict[params["problem"]["type"]](dom,farm,fs,bc)#,opt=opt)

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
    if params.get("optimization",{}):
        opt=windse.Optimizer(solver)
        if params["optimization"].get("taylor_test",False):
            opt.TaylorTest()

        if params["optimization"].get("optimize",True):
            opt.Optimize()

    tock = time.time()

    print("Run Complete: {:1.2f} s".format(tock-tick))

def main():
    actions = {"run": run_action}
    actions[get_action()]()

if __name__ == "__main__":
    main()
