from dolfin import *
from windse import windse_parameters
import numpy as np

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *
    from windse.blocks import blockify, MpiEvalBlock

def mpi_eval(u, x, comm=MPI.comm_world):
    # dim = u.function_space().dim()
    # mesh = u.function_space().mesh

    # # only evaluate on the processor that owns the point
    # val = np.nan
    # point = Point(x)
    # if mesh.bounding_box_tree().compute_first_entity_collision(point) < mesh.num_cells():
    #     val = u(point)
    #     rank = comm.rank








    u.set_allow_extrapolation(False) # use to make sure we don't evaluate outside of a processor

    try:
        ux = u(np.array(x))
    except:
        ux = np.nan

    u_list = comm.allgather(ux)

    for val in u_list:
        if not np.isnan(val):
            ux = val
            break

    return Constant(ux)




# blockify functions as needed
if windse_parameters.dolfin_adjoint:
    block_kwargs = {
        "base_eval": mpi_eval
    }
    mpi_eval = blockify(mpi_eval,MpiEvalBlock,block_kwargs=block_kwargs)
