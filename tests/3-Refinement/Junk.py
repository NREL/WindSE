import windse

### Alias Parameters ###
params = windse.windse_parameters

### Set General Parameters ###
params["general"]["name"] = "Refinement_Test"

### Set Box Parameters ###
params["domain"]["type"]    = "box"
params["domain"]["x_range"] = [-1000.0, 1000.0]
params["domain"]["y_range"] = [-1000.0, 1000.0]
params["domain"]["z_range"] = [    0.0,  500.0]
params["domain"]["nx"]      = 48
params["domain"]["ny"]      = 48
params["domain"]["nz"]      = 12
# params["solver"]["wind_range"] = [0.,3.14]
params["solver"]["wind_range"] = [3.14/4.,3.14]

### Set Wind Farm Parameters ###
params["wind_farm"]["type"]      = "grid"
params["wind_farm"]["grid_rows"] = 3
params["wind_farm"]["grid_cols"] = 3
params["wind_farm"]["ex_x"]      = [-350.0,350.0]
params["wind_farm"]["ex_y"]      = [-350.0,350.0]
params["wind_farm"]["HH"]        = 75.0
params["wind_farm"]["RD"]        = 50.0
params["wind_farm"]["thickness"] = 20.0
params["wind_farm"]["yaw"]       = 0.0
params["wind_farm"]["axial"]     = 0.33

### Set Refinement Parameters ###
# params["refine"]["farm_num"]       = 1
# params["refine"]["farm_type"]      = "cylinder"
# params["refine"]["farm_factor"]    = 1.25
# params["refine"]["turbine_num"]    = 1
# params["refine"]["turbine_type"]   = "simple"
# params["refine"]["turbine_factor"] = 2
params["refine"]["refine_custom"]  = [
    [ "cylinder", [ [0,0,0], 750, 150 ]               ],
    [ "box",      [ [[-500,500],[-500,500],[0,150]] ] ],
    [ "simple",   [ 100 ]                             ],
    [ "tear",     [ 50, 0.7853 ]                      ]
]

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
