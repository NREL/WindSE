import os
import __main__

### Get the name of program importing this package ###
if hasattr(__main__,"__file__"):
    main_file = os.path.basename(__main__.__file__)
else:
    main_file = "ipython"

### This checks if we are just doing documentation ###
if not main_file in ["sphinx-build", "__main__.py"]:
    import numpy as np
    import time
    import os.path as osp
    import argparse
    import sys
    from pyadjoint import get_working_tape
    import cProfile
    import dolfin

# import os
# os.environ['OMP_NUM_THREADS'] = '1'


ALL_ACTIONS = ("run", "mesh")
help_msg = """
Available commands:

    run      run windse with a specified params file
    mesh     export the mesh generated from a param file

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

    ### Clean up other module references ###
    mods_to_remove = []
    for k in sys.modules.keys():
        # if ("windse" in k):
        if ("windse" in k or "dolfin_adjoint" in k or "fenics_adjoint" in k):
            mods_to_remove.append(k)

    for i in range(len(mods_to_remove)):
        del sys.modules[mods_to_remove[i]]

    ### Clean tape if available ###
    tape = get_working_tape()
    if tape is not None:
        tape.clear_tape()

    ### Import fresh version of windse ###
    import windse
    try:
        from .driver_functions import SetupSimulation
    except:
        from driver_functions import SetupSimulation

    ### Setup everything ###
    params, problem, solver = SetupSimulation(params_loc)

    ### run the solver ###
    solver.Solve()

    problem.farm.finalize_farm()

    ### Perform Optimization ###
    if params.performing_opt_calc:
        opt=windse.Optimizer(solver)
        if params["optimization"]["gradient"] or params["general"]["debug_mode"]:
            opt.Gradient()
        if params["optimization"]["taylor_test"]:
            opt.TaylorTest()
        if params["optimization"]["optimize"]:
            opt.Optimize()


    tock = time.time()
    runtime = tock-tick
    if params.rank == 0:
        print("Run Complete: {:1.2f} s".format(runtime))

    # Begin postprocessing routines
    data_to_write = params["postprocessing"]["write_floris_input"]

    if data_to_write is not None:
        windse.write_to_floris(data_to_write, solver)

    params.comm.Barrier()


    dolfin.list_timings(dolfin.TimingClear.clear, [dolfin.TimingType.wall])

    return runtime

def mesh_action(params_loc=None):
    import windse
    from .driver_functions import Initialize, BuildDomain

    params = Initialize()
    dom, farm = BuildDomain(params)

    dom.ExportMesh()

    print("Mesh Exported to: "+dom.params.folder+"mesh/exported_mesh/")

def test_demo(demo):
    try:
        runtime = run_action(params_loc=demo)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        return (False,e,exc_traceback)
    else:
        return (True,runtime,None)


def main():
    actions = {"run":  run_action,
               "mesh": mesh_action}
    actions[get_action()]()

    # list_timings(TimingClear.clear, [TimingType.wall])

    # first_arg = get_action()

    # if first_arg == 'run': 
    #     cprof_action = 'run_action()'
    # elif first_arg == 'mesh': 
    #     cprof_action = 'mesh_action()'

    # cProfile.run(cprof_action, 'test')

if __name__ == "__main__":
    main()
    