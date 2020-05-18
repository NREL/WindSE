import numpy as np
import windse
from windse_driver.new_driver import run_driver, run_model, initialize_analysis, setup_problem, solve_problem

# Get the various objects associated with the problem from the driver
# params, problem, solver = run_model(params_loc="2d_wind_farm_induction_factor_opt.yaml")
params = initialize_analysis(params_loc="2d_wind_farm_induction_factor_opt.yaml", comm=None)
problem = setup_problem(params)


# Compute the functional
# print('functional value = ', solver.J)
print('axial values = ', problem.farm.a)

# Lets mess around with
new_a = np.array([0.28224965, 0.30361011, 0.30368474, 0.33775237, 0.26041658, 0.26898364, 0.27798243, 0.32165405, 0.25050653])
perturb_arr = np.array([0.03224965, 0.05361011, 0.05368474, 0.08775237, 0.01041658, 0.01898364, 0.02798243, 0.07165405, 0.00050653])

# print('\naxial values = ', problem.farm.a)
problem.farm.a[:] = new_a # += perturb_arr
problem.params.full_farm.UpdateConstants()
problem.ComputeFunctional()
solver = solve_problem(params, problem)
# solver.J = 0.0

print('functional value = ', solver.J)


# Create the optimization objects
opt = windse.Optimizer(solver)

# Get the Gradient
grad_val = opt.Gradient()
print(repr(grad_val))
