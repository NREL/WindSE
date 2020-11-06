import numpy as np
import time
import os.path as osp
import argparse
import sys
import cProfile
from dolfin import MPI

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
    pr = cProfile.Profile()
    pr.enable()

    ### Clean up other module references ###
    mods_to_remove = []
    for k in sys.modules.keys():
        if ("windse" in k):
        # if ("windse" in k or "fenics" in k or "pyadjoint" in k or "dolfin" in k):
            mods_to_remove.append(k)
    for i in range(len(mods_to_remove)):
        del sys.modules[mods_to_remove[i]]

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

    ### Perform Optimization ###
    if params.dolfin_adjoint:
        opt=windse.Optimizer(solver)
        if params["optimization"]["gradient"]:
            opt.Gradient()
        if params["optimization"]["taylor_test"]:
            opt.TaylorTest()
        if params["optimization"]["optimize"]:
            opt.Optimize()

    tock = time.time()
    runtime = tock-tick
    print("Run Complete: {:1.2f} s".format(runtime))


    comm = MPI.comm_world
    rank = comm.Get_rank()
    num_procs = comm.Get_size()

    # Dump the profiling logs for each rank to a file
    pr.dump_stats('%sprofiling/profiling_%d.prof' % (params.folder, rank))
    pr.disable()

    return runtime

def mesh_action():
    raise NotImplementedError

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
    