import warnings
import numpy as np
from pystatreduce.quantity_of_interest import QuantityOfInterest
import windse
from windse_driver.new_driver import run_driver, run_model, initialize_analysis, setup_problem, solve_problem

# Plotting specific imports - Temporary
from mpl_toolkits import mplot3d
import matplotlib
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=150)

class WindFarm(QuantityOfInterest):
    def __init__(self, systemsize, param_file, initial_solve=True, mpi_comm=None):
        QuantityOfInterest.__init__(self, systemsize)

        # Construct the WindSE problem
        self.initial_solve = initial_solve
        if self.initial_solve:
            self.windse_params, self.windse_problem, self.windse_solver = run_model(params_loc=param_file, comm=mpi_comm)
        else:
            warnings.warn("The solver object has not been created. This QoI is a compatibility object cannot be used for function evaluation")
            self.windse_params = initialize_analysis(params_loc=param_file, comm=None)
            self.windse_problem = setup_problem(self.windse_params)

        assert self.systemsize == self.windse_problem.farm.numturbs

        if hasattr(self, 'windse_solver'):
            self.windse_opt = windse.Optimizer(self.windse_solver)
        else:
            pass

        # # Create the WindSE optimization object
        # self.windse_opt = windse.Optimizer(self.windse_solver)


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

    def create_solver_object(self):
        if hasattr(self, 'windse_solver'):
            pass
        else:
            self.windse_solver = solve_problem(self.windse_params, self.windse_problem)
            self.windse_opt = windse.Optimizer(self.windse_solver) # You need this to compute gradients


if __name__ == '__main__':

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    params_file = "./3x3/2d_wind_farm_induction_factor_opt.yaml"
    # params_file = "./pruf/2d_wind_farm_PRUF.yaml"
    n_turbines = 9
    windfarm = WindFarm(n_turbines, params_file, mpi_comm=comm)

    # windfarm.plot_eigenmodes(np.ones(n_turbines))

    rv_val = 0.25*np.ones(n_turbines)
    fval = windfarm.eval_QoI(rv_val, np.zeros(n_turbines))

    grad_val = windfarm.eval_QoIGradient(rv_val, np.zeros(n_turbines))

    # Print it
    print('fval = ', fval)
    print('grad_val = \n', grad_val)
