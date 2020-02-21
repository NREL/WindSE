import numpy as np
import windse
from windse_driver.new_driver import run_driver, run_model
import time

from mpi4py import MPI


use_serial = False
if use_serial:
    start_time = time.time()
    # Get the various objects associated with the problem from the driver
    params, problem, solver = run_model(params_loc="./2d_wind_farm_PRUF.yaml")
    print("problem attributes: ", problem.__dict__.keys())
    print('params attributes: ', params.__dict__.keys())

    # # print the coordinates
    # print('problem.dom.x: \n', params.full_farm.x)
    # print('problem.dom.y: \n', params.full_farm.y)


    # Compute the functional
    print('functional value = ', solver.J)
    print('time_elapsed = ', time.time()-start_time)
    # print('\naxial values = ', problem.farm.a)
    # print(problem.full_farm.axial)
else:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        params, problem, solver = run_model(params_loc="./2d_wind_farm_PRUF.yaml")
        # Compute the functional
        print('functional value = ', solver.J)
        print('time_elapsed = ', time.time()-start_time)
    else:
        pass
