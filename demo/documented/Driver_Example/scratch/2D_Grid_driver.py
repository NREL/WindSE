import windse
import pprint

windse.initialize("params.yaml")

dom = windse.RectangleDomain()

farm = windse.GridWindFarm(dom)
farm.Plot(False)

fs = windse.TaylorHoodFunctionSpace(dom)

bc = windse.UniformInflow(dom,fs,farm)

problem = windse.TaylorHoodProblem(dom,farm,fs,bc)
# print('problem attributes = ', pprint.pprint(problem.__dict__))

solver = windse.SteadySolver(problem)
solver.Solve()

optimizer = windse.Optimizer(solver)
# gradient =
