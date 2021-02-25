import numpy as np
import time
import os.path as osp
import argparse
import sys
from pyadjoint import get_working_tape

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
    from .driver_functions import SetupSimulation

    ### Setup everything ###
    params, problem, solver = SetupSimulation(params_loc)

    ### run the solver ###
    solver.Solve()

    ### Perform Optimization ###
    if params.dolfin_adjoint:
        opt=windse.Optimizer(solver)
        if params["optimization"]["gradient"] or params["general"]["debug_mode"]:
            opt.Gradient()
        if params["optimization"]["taylor_test"]:
            opt.TaylorTest()
        if params["optimization"]["optimize"]:
            opt.Optimize()

    tock = time.time()
    runtime = tock-tick
    print("Run Complete: {:1.2f} s".format(runtime))

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

if __name__ == "__main__":
    main()
