import windse

windse.initialize("params.yaml")

dom = windse.RectangleDomain()

farm = windse.GridWindFarm(dom)
farm.Plot(False)

fs = windse.TaylorHoodFunctionSpace(dom)

bc = windse.UniformInflow(dom,fs,farm)

problem = windse.TaylorHoodProblem(dom,farm,fs,bc)
solver = windse.SteadySolver(problem)
solver.Solve()
