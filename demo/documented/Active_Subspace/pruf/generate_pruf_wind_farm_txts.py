# This file creates *.txt files with wind farm information that can then be
# loaded through a *.yaml file to recreate wind farms. In this file this is done
# to create samples for the active subspace method implementation where samples
# from a joint distribution are used to to create these text files.

import numpy as np
import chaospy as cp
import windse
from windse_driver.new_driver import initialize_analysis, setup_problem

params_loc = './2d_wind_farm_PRUF.yaml'
params = initialize_analysis(params_loc=params_loc, comm=None)
problem = setup_problem(params)

# Setup the joint distribution
n_turbines = 137
use_uniform = True
if use_uniform:
    print("Using uniform distribution")
    lb = 0.1
    ub = 0.6
    dist = cp.Uniform(lower=lb, upper=ub)
    jdist = cp.Iid(dist, n_turbines)
else:
    print("Using normal distribution")
    mu = 0.33
    std_dev = 0.075
    dist = cp.Normal(mu=mu, sigma=std_dev)
    jdist = cp.Iid(dist, n_turbines)

# Get the samples
n_samples = 1000
samples = jdist.sample(n_samples)

# finally update the induction factor and save txts
for i in range(0, n_samples):
    problem.farm.a[:] = samples[:,i]
    problem.params.full_farm.UpdateConstants()
    problem.params.full_farm.SaveWindFarm(val=i, filename="windfarm_pruf_uniform")
