import sys
import numpy as np
import time
import os.path as osp
import argparse
import sys

# import os
# os.environ['OMP_NUM_THREADS'] = '1'


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
def run_action(params_loc=None):
    tick = time.time()

    ### unload windse if previously loaded ###
    if "windse" in sys.modules.keys():
        del sys.modules["windse"]
        mods_to_remove = []
        for k in sys.modules.keys():
            if "windse." in k:
                mods_to_remove.append(k)
        for i in range(len(mods_to_remove)):
            del sys.modules[mods_to_remove[i]]

    import windse
    
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
        if params["domain"].get("scaled",False):
            Sx = 1.0e-3
        else:
            Sx = 1.0
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
            dom.WarpSplit(warp_height*Sx,warp_percent)


        if refine_custom is not None:
            for refine_data in refine_custom:
                if refine_data[1] == "full":
                    dom.Refine(refine_data[0])
                else:
                    if refine_data[1] == "custom":
                        region = np.array(refine_data[2])*Sx
                    else:
                        region = farm.CalculateFarmRegion(refine_data[1],length=refine_data[2]*Sx)
                    dom.Refine(refine_data[0],region=region,region_type=refine_data[1])
        
        if farm_num > 0:
            region = farm.CalculateFarmRegion(farm_type,farm_factor,length=farm_radius*Sx)
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
               "log":windse.LogLayerInflow,
               "turbsim":windse.TurbSimInflow}
    bc = bc_dict[params["boundary_condition"]["vel_profile"]](dom,fs,farm)

    ### Generate the problem ###
    prob_dict = {"stabilized":windse.StabilizedProblem,
                 "taylor_hood":windse.TaylorHoodProblem,
                 "unsteady":windse.UnsteadyProblem}
    problem = prob_dict[params["problem"]["type"]](dom,farm,fs,bc)#,opt=opt)

    ### Solve ###
    solve_dict = {"steady":windse.SteadySolver,
                  "unsteady":windse.UnsteadySolver,
                  "multiangle":windse.MultiAngleSolver,
                  "importedvelocity":windse.TimeSeriesSolver}
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

        if params["optimization"].get("optimize",False):
            opt.Optimize()

        if params["optimization"].get("gradient",True):
            opt.Gradient()

    tock = time.time()
    runtime = tock-tick
    print("Run Complete: {:1.2f} s".format(runtime))

    return runtime


def test_demo(demo):
    try:
        runtime = run_action(params_loc=demo)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        return (False,e,exc_traceback)
    else:
        return (True,runtime,None)


def main():
    actions = {"run":  run_action}
    actions[get_action()]()

if __name__ == "__main__":
    main()
