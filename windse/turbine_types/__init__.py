from dolfin import *
from windse import windse_parameters

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *

### Import the turbine types
from .GenericTurbine  import GenericTurbine
from .Actuator2D      import Actuator2D
from .ActuatorDisk    import ActuatorDisk
from .ActuatorLine    import ActuatorLine

### Create the turbine dictionary ###
turbine_dict = {
    "disk2D": Actuator2D,
    "disk":   ActuatorDisk,
    "line":   ActuatorLine,
}