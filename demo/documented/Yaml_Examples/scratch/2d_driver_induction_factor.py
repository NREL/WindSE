import windse
from windse_driver.new_driver import run_driver, run_model

# Get the various objects associated with the problem from the driver
params, problem, solver = run_model(params_loc="2d_wind_farm_induction_factor_opt.yaml")

# Compute the functional
print('functional value = ', solver.J)
print(problem.__dict__.keys())
print(params.__dict__.keys())
print('\naxial values = ', problem.farm.a)
# print(problem.full_farm.axial)
print('type of problem:', type(problem))
print('\n', dir(problem))
# Lets do the optimization
# Create the optimization objects
opt = windse.Optimizer(solver)

# Get the Gradient
grad_val = opt.Gradient()
print(type(grad_val))
