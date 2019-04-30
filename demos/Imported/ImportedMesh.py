import windse

### Initialize WindSE ###
windse.initialize("params.yaml")

### Generate Domain ###
dom = windse.ImportedDomain()

### Generate Wind Farm ###
farm = windse.GridWindFarm(dom)
farm.Plot()

### Refine the Domain ###
# region = [[farm.ex_x[0],dom.x_range[1]],farm.ex_y,farm.ex_z]
# region = [[-1000,2500],[-1000,1000],[0,216]]
# dom.Refine(1,region=region)
dom.Save()

### Function Space ###
fs = windse.LinearFunctionSpace(dom)
print(fs.W.dim())

### Setup Boundary Conditions ###
bc = windse.PowerInflow(dom,fs)

### Generate the problem ###
problem = windse.StabilizedProblem(dom,farm,fs,bc)

### Solve ###
solver = windse.SteadySolver(problem)
solver.Solve()

### Output Results ###
solver.Save()