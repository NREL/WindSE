import windse_driver.driver_functions as df

#####################################


### Load the parameters from file ###
# windse.initialize("params.yaml")
# params = windse.windse_parameters





### Manually define the parameters ###
### Create an default parameter set ###
params = df.BlankParameters()

### Edit the general settings ###
params["general"]["name"]        = "2D_driver"
params["general"]["output"]      = ["mesh","initial_guess","turbine_force","solution"]
params["general"]["output_type"] = "xdmf"

### Edit the windfarm settings ###
params["wind_farm"]["type"]      = "grid"
params["wind_farm"]["grid_rows"] = 6
params["wind_farm"]["grid_cols"] = 6
params["wind_farm"]["ex_x"]      = [-1800,1800]
params["wind_farm"]["ex_y"]      = [-1800,1800]

### Edit the turbine settings
params["turbines"]["type"]      = "2D_disk"
params["turbines"]["HH"]        = 90
params["turbines"]["RD"]        = 126
params["turbines"]["thickness"] = 10
params["turbines"]["yaw"]       = 0
params["turbines"]["axial"]     = 0.33

### Edit the domain settings ###
params["domain"]["type"]    = "rectangle"
params["domain"]["x_range"] = [-2500, 2500]
params["domain"]["y_range"] = [-2500, 2500]
params["domain"]["nx"]      = 50
params["domain"]["ny"]      = 50

### Edit the boundary condition settings ###
params["boundary_conditions"]["vel_profile"] = "uniform"

### Edit the problem settings ###
params["function_space"]["type"] = "taylor_hood"
params["problem"]["type"]        = "taylor_hood"
params["solver"]["type"]         = "steady"

### Initialize the parameters object ###
params = df.Initialize(params)

### Build the domain and wind-farm ###
dom, farm = df.BuildDomain(params)
dom.Save()

### Compile the problem and Build the solver ###
problem = df.BuildProblem(params,dom,farm)
solver = df.BuildSolver(params,problem)

### Solve and save ###
solver.Solve()
solver.Save()
