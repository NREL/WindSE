from dolfin import *
from windse import windse_parameters
import numpy as np
import time

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *
    from windse.blocks import blockify, MpiEvalBlock, UflEvalBlock
    from pyadjoint.overloaded_type import create_overloaded_object

def meteor_to_math(angle):
    '''
    this function converts an angle from the meteorological convention to the math convention.
    0   -> 270
    90  -> 180
    180 -> 90
    270 -> 0
    Additionally the input is in degrees and the output is in radians. 

    '''
    return np.radians(270.0-angle)

def math_to_meteor(angle):
    '''
    this function converts an angle from the math convention to the meteorological convention.
    0   -> 270
    90  -> 180
    180 -> 90
    270 -> 0
    Additionally the input is in radians and the output is in degrees. 

    '''
    return 270.0-np.degrees(angle)

def ufl_eval(form, print_statement=None):
    '''
    This function converts complex ufl forms to floats
    '''
    # mesh = UnitCubeMesh(2,2,2) # TODO: this might be bad in parallel
    # dx_lame = Measure("dx",mesh)
    # out = assemble(form*dx_lame)
    # return out
    # if print_statement is not None:
    #     print(f'Standard: {print_statement}, output: {float(Constant(form)):1.20f}')
    return Constant(form, name = "ufl_eval")

if windse_parameters.dolfin_adjoint:
    block_kwargs = {
        "base_eval": ufl_eval
    }
    ufl_eval = blockify(ufl_eval,UflEvalBlock,block_kwargs=block_kwargs)


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
    


    assert len(x) == u.geometric_dimension()



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
    # if mesh.bounding_box_tree().compute_first_entity_collision(point) < mesh.num_cells():
    #     ux = np.array(u(point))

    try:
        ux = np.array(u(point))
        # print(np.array(x,dtype=float))
        # print(f"proc: {rank} owns {ufl_eval(x).values()}")
    except:
        pass
        # print("fail: "+repr(rank))


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

    # if value_rank == 0:
    #     out = Constant(ux)
    # elif value_rank == 1:
    #     out = numpy_adjoint.array.ndarray(ux)
    # print(type(ux))
    if windse_parameters.dolfin_adjoint:
        out =  create_overloaded_object(ux)
    else:
        out = ux
        
    # out = Constant(ux)
    # print(type(out))

    # print(f"x: {np.array(x,dtype=float)}")
    # print(f"x: {np.array(x,dtype=float)}")
    # print(f"u: {np.array(out,dtype=float)}")
    # print(type(out))
    return out

if windse_parameters.dolfin_adjoint:
    block_kwargs = {
        "base_eval": mpi_eval
    }
    mpi_eval = blockify(mpi_eval,MpiEvalBlock,block_kwargs=block_kwargs)



def test_dolfin_adjoint(control_list,form):
    tick = time.time()
    if not isinstance(form, AdjFloat):
        if form.ufl_domain() is None:
            mesh = UnitCubeMesh(8,8,8)
            dx = Measure("dx",mesh)
            J = assemble(form*dx)
        else:
            dx = Measure("dx",form.ufl_domain())
            J = assemble(form*dx)
    else:
        J = form

    h = []
    controls = []
    init_vals = []
    for c in control_list:
        h.append(Constant(0.1*float(c)))
        controls.append(Control(c))
        init_vals.append(Constant(float(c)))

    Jhat = ReducedFunctional(J,controls)

    conv_rate = taylor_test(Jhat, init_vals, h)

    der = Jhat.derivative()

    for d in der:
        print(d.values())
    tock = time.time()
    print(f"Total time: {tock-tick:1.2f} s")






