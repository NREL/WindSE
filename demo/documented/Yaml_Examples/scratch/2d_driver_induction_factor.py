import numpy as np
import windse
from windse_driver.new_driver import run_driver, run_model

# Get the various objects associated with the problem from the driver
params, problem, solver = run_model(params_loc="2d_wind_farm_induction_factor_opt.yaml")
print("problem attributes: ", problem.__dict__.keys())
print('params attributes: ', params.__dict__.keys())

# # print the coordinates
# print('problem.dom.x: \n', params.full_farm.x)
# print('problem.dom.y: \n', params.full_farm.y)


# Compute the functional
print('functional value = ', solver.J)
# print('\naxial values = ', problem.farm.a)
# print(problem.full_farm.axial)


# Lets mess around with
perturb_arr = np.array([0.03224965, 0.05361011, 0.05368474, 0.08775237, 0.01041658, 0.01898364, 0.02798243, 0.07165405, 0.00050653])
problem.farm.a[:] += perturb_arr
print('\naxial values = ', problem.farm.a)
problem.params.full_farm.UpdateConstants()
# solver.J = 0.0
# problem.ComputeFunctional()
# solver.J = 0.0
solver.Solve()

print('functional value = ', solver.J)


# Lets do the optimization
# Create the optimization objects
opt = windse.Optimizer(solver)

# Get the Gradient
grad_val = opt.Gradient()
print(type(grad_val))
