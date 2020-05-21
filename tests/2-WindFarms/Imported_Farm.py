import windse
import numpy as np

### Alias Parameters ###
params = windse.windse_parameters

### Set General Parameters ###
params["general"]["name"] = "Grid_Farm_Test"

### Set Box Parameters ###
params["domain"]["type"]    = "box"
params["domain"]["x_range"] = [-1000.0, 1000.0]
params["domain"]["y_range"] = [-1000.0, 1000.0]
params["domain"]["z_range"] = [    0.0,  500.0]
params["domain"]["nx"]      = 24
params["domain"]["ny"]      = 24
params["domain"]["nz"]      = 6

### Set Wind Farm Parameters ###
params["wind_farm"]["type"]      = "grid"
params["wind_farm"]["grid_rows"] = 4
params["wind_farm"]["grid_cols"] = 3
params["wind_farm"]["ex_x"]      = [-300.0,300.0]
params["wind_farm"]["ex_y"]      = [-300.0,300.0]
params["wind_farm"]["HH"]        = 90.0
params["wind_farm"]["RD"]        = 126.0
params["wind_farm"]["thickness"] = 20.0
params["wind_farm"]["yaw"]       = 0.0
params["wind_farm"]["axial"]     = 0.33

### Initialize Parameters using those set above ###
windse.initialize(None)

### Create the Domain Object ###
dom = windse.BoxDomain()

### Create Wind Farm Object ###
farm = windse.GridWindFarm(dom)

### Refine Domain ###
windse.RefineMesh(dom,farm)

### Check if the object is as expected ###
dom.Save()
farm.Plot(False)

### Set True Values ###
x = [-237.,    0.,  237., -237.,    0.,  237., -237.,    0.,  237., -237.,    0.,  237.]
y = [-237., -237., -237.,  -79.,  -79.,  -79.,   79.,   79.,   79.,  237.,  237.,  237.]
z = [90., 90., 90., 90., 90., 90., 90., 90., 90., 90., 90., 90.]


### Check True Values ###
if np.linalg.norm(x-farm.x) > 1e-10:
    print("Expected x-values: " + repr(x))
    print("Actual x-values:   " + repr(list(farm.x)))
    raise ValueError("Grid Farm constructed with unexpected x locations")
if np.linalg.norm(y-farm.y) > 1e-10:
    print("Expected x-values: " + repr(y))
    print("Actual x-values:   " + repr(list(farm.y)))
    raise ValueError("Grid Farm constructed with unexpected y locations")
if np.linalg.norm(z-farm.z) > 1e-10:
    print("Expected x-values: " + repr(z))
    print("Actual x-values:   " + repr(list(farm.z)))
    raise ValueError("Grid Farm constructed with unexpected z locations")
