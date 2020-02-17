import numpy as np
from pystatreduce.quantity_of_interest import QuantityOfInterest
import windse
from windse_driver.new_driver import run_driver, run_model

# Plotting specific imports - Temporary
from mpl_toolkits import mplot3d
import matplotlib
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=150)

class WindFarm(QuantityOfInterest):
    def __init__(self, systemsize, param_file):
        QuantityOfInterest.__init__(self, systemsize)

        # Construct the WindSE problem
        self.windse_params, self.windse_problem, self.windse_solver = run_model(params_loc=param_file)
        assert self.systemsize == self.windse_problem.farm.numturbs

        # Create the WindSE optimization object
        self.windse_opt = windse.Optimizer(self.windse_solver)


    def eval_QoI(self, mu, xi):
        rv = mu + xi
        self.windse_problem.farm.a[:] = rv
        self.windse_problem.params.full_farm.UpdateConstants()
        self.windse_problem.ComputeFunctional()
        self.windse_solver.Solve()
        return self.windse_solver.J

    def eval_QoIGradient(self, mu, xi):
        rv = mu + xi
        self.windse_problem.farm.a[:] = rv
        self.windse_problem.params.full_farm.UpdateConstants()
        self.windse_problem.ComputeFunctional()
        self.windse_solver.Solve()
        return self.windse_opt.Gradient()

    def plot_eigenmodes(self, eigenvec, is_2D=True):
        """
        Plotting script for ploting eigenvectors on the windfarm. This is a
        temporary home for this function as it doesnt really belong in here.
        """
        # Get the wind farm bounds
        x, y, z = self.windse_problem.farm.GetLocations()
        fname = 'eigenvector_overlay.pdf'
        fig = plt.figure('eigenvectors')
        ax = plt.axes()
        s = ax.scatter(x, y, c=eigenvec, cmap="coolwarm", edgecolors=(0,0,0,1))
        plt.tight_layout()
        plt.show()

if __name__ == '__main__':

    params_file = "./2d_wind_farm_induction_factor_opt.yaml"
    n_turbines = 9
    windfarm = WindFarm(n_turbines, params_file)

    windfarm.plot_eigenmodes(np.ones(n_turbines))

    """
    rv_val = 0.25*np.ones(n_turbines)
    fval = windfarm.eval_QoI(rv_val, np.zeros(n_turbines))

    grad_val = windfarm.eval_QoIGradient(rv_val, np.zeros(n_turbines))

    # Print it
    print('fval = ', fval)
    print('grad_val = \n', grad_val)
    """
