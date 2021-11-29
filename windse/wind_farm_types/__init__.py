### Import dolfin
from dolfin import *
from windse import windse_parameters
from windse.turbine_types import turbine_dict

### Check if we need dolfin_adjoint ###
if windse_parameters.dolfin_adjoint:
    from dolfin_adjoint import *
    
### import the wind farm types
from .GenericWindFarm  import GenericWindFarm
from .GridWindFarm     import GridWindFarm
from .RandomWindFarm   import RandomWindFarm
from .ImportedWindFarm import ImportedWindFarm
from .EmptyWindFarm    import EmptyWindFarm

