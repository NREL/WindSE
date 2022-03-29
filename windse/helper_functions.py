from dolfin import *
from windse import windse_parameters
import numpy as np

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *
    from windse.blocks import blockify, MpiEvalBlock

def mpi_eval(u, x, comm=MPI.comm_world):

    '''

    This function returns the value of a function, ``u``, at a point, ``x``. 
    It is written to carry out this operation in parallel where not all
    processes will have ownership of the point ``x``. In this case, ownership 
    is established with ``bounding_box_tree().compute_first_entity_collision(point)``
    where the resulting interpolated value is then distributed to all processes 
    immediately afterward.
    
    Args:
        u (dolfin function):
            The dolfin function to be sampled. 
        x ([1x2], [1x3] list):
            The point at which the function should be sampled, must be within the domain of ``u``. 

    Returns:
        Constant(ux):
            The interpolated value of function ``u`` at point ``x`` as a dolfin Constant. 

    '''

    rank = comm.Get_rank()
    num_procs = comm.Get_size()

    value_rank = u.value_rank()
    mesh = u.function_space().mesh()
    u.set_allow_extrapolation(False)

    point = Point(x)

    if value_rank == 0:
        # For a scalar function, create a single placeholder
        ux = np.array(np.nan)
        nn = 1

    elif value_rank == 1:
        # For a vector function, create a placeholder with ndim entries
        ux = np.array([np.nan for k in range(u.geometric_dimension())])
        nn = len(ux)

    ux_global = np.zeros(nn*num_procs)

    # Check if this rank owns the point and, if so, evaluate the function there 
    if mesh.bounding_box_tree().compute_first_entity_collision(point) < mesh.num_cells():
        ux = np.array(u(point))


    # implementation_method = 1
    implementation_method = 2 # 2 is recommended, type(ux) is variable with method 1

    if implementation_method == 1:

        comm.Gather(ux, ux_global, root=0)

        if rank == 0:
            ux_global = ux_global.reshape(-1, nn)
            ux = np.nanmean(ux_global, axis=0)

        comm.Bcast(ux, root=0)

    elif implementation_method == 2:

        comm.Allgather(ux, ux_global)

        ux_global = ux_global.reshape(-1, nn)
        ux = np.nanmean(ux_global, axis=0)

    return Constant(ux)




# blockify functions as needed
if windse_parameters.dolfin_adjoint:
    block_kwargs = {
        "base_eval": mpi_eval
    }
    mpi_eval = blockify(mpi_eval,MpiEvalBlock,block_kwargs=block_kwargs)
