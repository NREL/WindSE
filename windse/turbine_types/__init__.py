from dolfin import *
from windse import windse_parameters

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *

### Import the turbine types
from .GenericTurbine     import GenericTurbine
from .ActuatorDisk       import ActuatorDisk
from .ActuatorDisk2D     import ActuatorDisk2D
from .ActuatorDiskExpr     import ActuatorDiskExpr
from .ActuatorDiskNumpy  import ActuatorDiskNumpy
from .ActuatorLine       import ActuatorLine
from .ActuatorLineDolfin import ActuatorLineDolfin
from .DisabledTurbine    import DisabledTurbine

### Create the turbine dictionary ###
turbine_dict = {
    "disk":         ActuatorDisk,
    "disks":        ActuatorDisk,
    "2D_disk":      ActuatorDisk2D,
    "2D_disks":     ActuatorDisk2D,
    "expr_disk": ActuatorDiskExpr,
    "expr_disks": ActuatorDiskExpr,
    "numpy_disk":   ActuatorDiskNumpy,
    "numpy_disks":  ActuatorDiskNumpy,
    "line":         ActuatorLine,
    "lines":        ActuatorLine,
    "dolfin_line":  ActuatorLineDolfin,
    "dolfin_lines": ActuatorLineDolfin,
    "disabled":       DisabledTurbine,
}