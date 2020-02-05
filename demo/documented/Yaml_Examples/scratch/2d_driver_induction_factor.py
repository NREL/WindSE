import windse
from windse_driver.new_driver import run_driver

# Get the various objects associated with the problem from the driver
problem, params, solver = run_model(params_loc="2d_wind_farm_induction_factor_opt.yaml")

# Compute the functional
print('functional value = ', solver.J)

# Lets do the optimization
# Create the optimization objects
opt = windse.Optimizer(solver)

# Get the Gradient
grad_val = opt.Gradient()
print(type(grad_val))
