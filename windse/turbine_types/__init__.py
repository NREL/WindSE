from dolfin import *
from windse import windse_parameters

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *

### Import the turbine types
from .GenericTurbine       import GenericTurbine
from .ActuatorDisk         import ActuatorDisk
from .ActuatorDisk2D       import ActuatorDisk2D
from .ActuatorDiskNumpy    import ActuatorDiskNumpy
from .ActuatorLine         import ActuatorLine

### Create the turbine dictionary ###
turbine_dict = {
    "disk":       ActuatorDisk,
    "2D_disk":    ActuatorDisk2D,
    "numpy_disk": ActuatorDiskNumpy,
    "line":       ActuatorLine,
}