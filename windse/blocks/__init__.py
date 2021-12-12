from dolfin import *
from windse import windse_parameters
from pyadjoint.tape import get_working_tape, annotate_tape, stop_annotating
from pyadjoint.overloaded_type import create_overloaded_object

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *

### Import the turbine types
from .ActuatorDiskForceBlock import ActuatorDiskForceBlock
from .ActuatorLineForceBlock import ActuatorLineForceBlock
from .BaseHeightBlock        import BaseHeightBlock
from .ControlUpdaterBlock    import ControlUpdaterBlock
from .MarkerBlock            import MarkerBlock

### A function for automatically creating blocks
def blockify(function_backend, function_block, block_kwargs={}):
    '''
    this method helps dolfin adjoint track custom function by linking said function to a Block object

    Args: 
        function_backend (method):                   the function we want to wrap.
        function_block   (:class:`pyadjoint.Block`): the block class specific to this function that help build derivative
    '''
    
    ### Create a wrapper around the function we are blockifying
    def blockified(*args, **kwargs):
        
        ### Check if we are annotating 
        annotate = annotate_tape(kwargs)

        ### Evaluate the original function
        with stop_annotating():
            output = function_backend(*args,**kwargs)
        output = create_overloaded_object(output)

        if annotate:        
            block = function_block(*args, **kwargs, **block_kwargs)
            block.add_output(output.create_block_variable())

            tape = get_working_tape()
            tape.add_block(block)

        return output

    return blockified