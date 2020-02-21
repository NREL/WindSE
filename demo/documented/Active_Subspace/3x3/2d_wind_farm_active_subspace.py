# This file contains a script for running active subspace on the windfarm using
# pyStatReduce. We are still calling chaospy directly
import numpy as np
import chaospy as cp
import windse
from windse_driver.new_driver import run_driver, run_model
from demo.documented.Yaml_Examples.scratch.wind_farm_qoi import WindFarm
from pystatreduce.active_subspace import ActiveSubspace
from pystatreduce.dimension_reduction import DimensionReduction
import time

np.printoptions(linewidth=150)

run_analysis = False
plot_eigenvec = True

n_turbines = 9
use_uniform = False
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

# Compute the mean and variance
print('mean = ', cp.E(jdist))
print('std_dev = ', cp.Std(jdist))

# Construct the QOI
params_file = "./2d_wind_farm_induction_factor_opt.yaml"
n_turbines = 9
windfarm = WindFarm(n_turbines, params_file)

if run_analysis:
    start_time = time.time()

    use_active_subspace = False
    if use_active_subspace:
        use_arnoldi_sampling = False
        # Create the active subspace object
        n_samples = 1000
        active_subspace = ActiveSubspace(windfarm, n_dominant_dimensions=9, n_monte_carlo_samples=n_samples)
        active_subspace.getDominantDirections(windfarm, jdist)
        eigenvals = active_subspace.iso_eigenvals
        eigenvecs = active_subspace.iso_eigenvecs

        # Fname for saving
        if use_uniform:
            fname = './eigenmodes_active_uniform_' + str(n_samples)
        else:
            fname = './eigenmodes_active_normal_' + str(n_samples)
    else:
        use_arnoldi_sampling = True
        n_arnoldi_sample = n_turbines+1
        dominant_space = DimensionReduction(n_arnoldi_sample=n_arnoldi_sample, # systemsize+1,
                                        exact_Hessian=False,
                                        sample_radius=1.e-6)
        dominant_space.getDominantDirections(windfarm, jdist, max_eigenmodes=10)
        eigenvals = dominant_space.iso_eigenvals
        eigenvecs = dominant_space.iso_eigenvecs

        # Fname for saving
        if use_uniform:
            fname = './eigenmodes_dominant_uniform'
        else:
            fname = './eigenmodes_dominant_normal'

    # Print the eigenvalues and eigenvectors
    print('eigenvals = \n', eigenvals)
    print('\neigenvecs = \n', eigenvecs)

    # Finally save an npz file
    np.savez(fname, eigenvals=eigenvals, eigenvecs=eigenvecs)

    end_time = time.time()
    time_elapsed = end_time - start_time
    print('time elapsed = ', time_elapsed)

if plot_eigenvec:
    # Read the eigenmodes
    fname = "eigenmodes_active_uniform_1000.npz"
    eigenmodes = np.load(fname)
    eigenvecs = eigenmodes['eigenvecs']
    windfarm.plot_eigenmodes(eigenvecs[:,0])
