# This file aggregates the gradients from the ./gradients directory into a 2d
# array that is accepted by the Active Subspace method.
import glob
import numpy as np
import chaospy as cp
from pystatreduce.active_subspace import ActiveSubspace, BifidelityActiveSubspace
import pystatreduce.utils as utils

# Plotting specific imports
from mpl_toolkits import mplot3d
import matplotlib
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=150)

# Boolean Options
read_case_npzs = False # REad in npz files. needed if you want to compute active subspace
compute_active_subspace = True # Bool to either read a file with eigenmodes in it
use_single_fidelity = False
compute_subspace_angles = False
plot_eigenvals = True # Plots the eigenvalues
plot_shadow_plot = False # Plots the shadow plot
get_trendline = False

n_turbines = 137
distribution_type = "uniform"

# Read the gradients
gradient_directory = "./gradients/"
n_samples = len(glob.glob1(gradient_directory, "*.npz"))
gradient_array = np.zeros([n_samples, n_turbines])
power_array = np.zeros(n_samples)
induction_factor_array = np.zeros([n_samples, n_turbines])

if read_case_npzs:
    for i in range(0, n_samples):
        fpath = gradient_directory + 'grad_2d_wind_farm_PRUF_' + distribution_type + '_' + str(i) + '.npz'
        f = np.load(fpath)
        gradient_array[i,:] = f['dJdinduction']
        power_array[i] = f['J']
        induction_factor_array[i,:] = f['a']

    # Uncomment below if you want to store the aggregated values into an npz file
    # np.savez('./lofi_results.npz', gradients=gradient_array, power_vals=power_array, induction_factors=induction_factor_array)


if compute_active_subspace:
    from demo.documented.Active_Subspace.wind_farm_qoi import WindFarm
    # Create the joint distribution
    if distribution_type is "uniform":
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

    # Create the QoI for Feeding it to the active subspace. This is not really
    # needed for computation, rather to be API compliant.
    nominal_yaml_fname = './2d_wind_farm_PRUF.yaml'
    pruf_qoi = WindFarm(n_turbines, nominal_yaml_fname, initial_solve=False)


    # Create the active subspace object and get the eigenmodes
    if use_single_fidelity:
        print('\n*** Using Single Fidelity ***\n')
        fpath = "./hifi_results.npz"
        f = np.load(fpath)
        gradient_array = f['gradients']
        power_array = f['power_vals']
        induction_factor_array = f['induction_factors']
        active_subspace = ActiveSubspace(pruf_qoi, n_dominant_dimensions=10,
                                         n_monte_carlo_samples=n_samples,
                                         read_gradient_samples=True,
                                         gradient_array=gradient_array)
        active_subspace.getDominantDirections(pruf_qoi, jdist)
        eigenmode_fname = './pruf_active_eigenmodes_unifi.npz'
    else: # Use bifidelity
        QoI_dict = {"high_fidelity" : pruf_qoi,
                    "low_fidelity" : pruf_qoi}
        # Read in the aggregated dictionaries that are generated
        lofi_fpath = "./lofi_results.npz"
        hifi_fpath = "./hifi_results.npz"
        split = [150, 1000] # Hifi vs lofi cases, delta is lofi - hifi
        n_active_subspace_samples = [150, 850]
        f1 = np.load(lofi_fpath)
        f2 = np.load(hifi_fpath)
        gradient_dict = {'high_fidelity' : -f2['gradients'][0:150,:],
                         'low_fidelity' : -f1['gradients']}
        active_subspace = BifidelityActiveSubspace(n_turbines, n_dominant_dimensions=1,
                                                   n_monte_carlo_samples=n_active_subspace_samples,
                                                   read_rv_samples=False, write_rv_samples=False,
                                                   read_gradient_samples=True,
                                                   gradient_dict=gradient_dict)
        active_subspace.getDominantDirections(QoI_dict, jdist)
        eigenmode_fname = './pruf_active_eigenmodes_bifi.npz'

    # Print the eigenmodes
    print('iso_eigenvals = \n', active_subspace.iso_eigenvals)
    print('iso_eigenvecs = \n', active_subspace.iso_eigenvecs)

    np.savez(eigenmode_fname, eigenvals=active_subspace.iso_eigenvals, eigenvecs=active_subspace.iso_eigenvecs)

    eigenvals = active_subspace.iso_eigenvals
    eigenvecs = active_subspace.iso_eigenvecs

    # Plot the most prominent eigenvector on the windfarm
    pruf_qoi.plot_eigenvector_on_farm(-eigenvecs[:,0], show=False)

if compute_subspace_angles:
    eigenmode_path1 = 'pruf_active_eigenmodes_unifi.npz'
    eigenmode_path2 = 'pruf_active_eigenmodes_bifi.npz'
    unifi_eigenmodes = np.load(eigenmode_path1)
    bifi_eigenmodes = np.load(eigenmode_path2)
    V1 = unifi_eigenmodes['eigenvecs'][:,0]
    V2 = bifi_eigenmodes['eigenvecs'][:,0]
    angles = utils.compute_subspace_angles(V1, V2)
    print('angles = \n', np.degrees(angles))

if plot_eigenvals:

    n_eigenmodes=10
    xrange = np.arange(n_eigenmodes)
    normalized_eigenvals = eigenvals / eigenvals[0]

    matplotlib.rcParams.update({'font.size': 18})
    fname = './pruf_active_eigenvalues_' + str(n_eigenmodes) + '.pdf'
    fig = plt.figure('eigenvalues', figsize=(15,6))
    ax = plt.axes()
    p1 = ax.plot(xrange, normalized_eigenvals[0:n_eigenmodes], color='k')
    if n_eigenmodes < 15:
        p2 = ax.plot(xrange, normalized_eigenvals[0:n_eigenmodes], 'o', color='k')
        ax.set_xticks(np.arange(n_eigenmodes))
    ax.set_yscale('log')
    ax.set_xlabel('index')
    ax.set_ylabel('normalized eigenvalues')
    plt.tight_layout()
    # plt.show()
    plt.savefig(fname)


if plot_shadow_plot:
    if not compute_active_subspace:
        if use_single_fidelity:
            eigenmode_path = 'pruf_active_eigenmodes_unifi.npz'
            f = np.load(eigenmode_path)
            eigenvecs = f['eigenvecs']
            fpath = "./hifi_results.npz"
        else: # Implement bifidelity
            eigenmode_path = 'pruf_active_eigenmodes_bifi.npz'
            bifi_eigenmodes = np.load(eigenmode_path)
            eigenvecs = bifi_eigenmodes['eigenvecs']
            fpath = "./lofi_results.npz" # Currently we will only consider the lo-fi results

        # Load the results for the different runs
        f = np.load(fpath)
        gradient_array = f['gradients']
        power_array = f['power_vals']
        induction_factor_array = f['induction_factors']

    eigenvector_index = 0
    x_arr = np.zeros(n_samples)
    for i in range(n_samples):
        x_arr[i] = np.dot(eigenvecs[:,0], induction_factor_array[i,:])

    # Normalize x_arr about its mean
    x_mean = np.mean(x_arr)
    x_arr -= x_mean

    # Finally plot the information
    if use_single_fidelity:
        fname = 'pruf_shadow_plot_unifi.pdf'
    else:
        fname = 'pruf_shadow_plot_bifi_coarse.pdf'
    fig = plt.figure('shadow_plot')
    ax = plt.axes()
    s = ax.scatter(x_arr, -power_array, edgecolors=(0, 0, 0, 1))
    ax.set_xlabel('W1*x')
    ax.set_ylabel('Power')
    plt.tight_layout()
    # plt.show()
    plt.savefig(fname)

if get_trendline:
    from scipy import stats

    if use_single_fidelity:
        eigenmode_path = 'pruf_active_eigenmodes_unifi.npz'
        f = np.load(eigenmode_path)
        eigenvecs = f['eigenvecs']
        fpath = "./hifi_results.npz"
    else: # Implement bifidelity
        eigenmode_path = 'pruf_active_eigenmodes_bifi.npz'
        bifi_eigenmodes = np.load(eigenmode_path)
        eigenvecs = bifi_eigenmodes['eigenvecs']
        fpath = "./hifi_results.npz" # Currently we will only consider the lo-fi results

    # Load the results for the different runs
    f = np.load(fpath)
    gradient_array = f['gradients']
    power_array = f['power_vals']
    induction_factor_array = f['induction_factors']

    eigenvector_index = 0
    x_arr = np.zeros(n_samples)
    for i in range(n_samples):
        x_arr[i] = np.dot(eigenvecs[:,0], induction_factor_array[i,:])

    # Normalize x_arr about its mean
    x_mean = np.mean(x_arr)
    x_arr -= x_mean

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_arr, -power_array)
    print('slope = ', slope)
    print('intercept = ', intercept)
    print('r_value = ', r_value)
    print('p_value = ', p_value)
    print('std_err = ', std_err)
