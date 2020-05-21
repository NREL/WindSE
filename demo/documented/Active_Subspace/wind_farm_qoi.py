import warnings
import numpy as np
from pystatreduce.quantity_of_interest import QuantityOfInterest
import windse
from windse_driver.driver_functions import SetupSimulation
# from windse_driver.new_driver import run_driver, run_model, initialize_analysis, setup_problem, solve_problem

# Plotting specific imports - Temporary
from mpl_toolkits import mplot3d
import matplotlib
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=150)

class WindFarm(QuantityOfInterest):
    def __init__(self, systemsize, param_file, initial_solve=True, mpi_comm=None):
        QuantityOfInterest.__init__(self, systemsize)

        # Get the necessary attributes
        self.windse_params, self.windse_problem, self.windse_solver = SetupSimulation(param_file)

        if initial_solve:
            self.windse_solver.Solve()

        if self.windse_params.dolfin_adjoint:
            self.windse_opt = windse.Optimizer(self.windse_solver)


        ##### OLD CODE #####
        # # Construct the WindSE problem
        # self.initial_solve = initial_solve
        # if self.initial_solve:
        #     self.windse_params, self.windse_problem, self.windse_solver = run_model(params_loc=param_file, comm=mpi_comm)
        # else:
        #     warnings.warn("The solver object has not been created. This QoI is a compatibility object cannot be used for function evaluation")
        #     self.windse_params = initialize_analysis(params_loc=param_file, comm=None)
        #     self.windse_problem = setup_problem(self.windse_params)
        # assert self.systemsize == self.windse_problem.farm.numturbs
        # if hasattr(self, 'windse_solver'):
        #     self.windse_opt = windse.Optimizer(self.windse_solver)
        # else:
        #     pass


    def eval_QoI(self, mu, xi):
        rv = mu + xi
        self.windse_problem.farm.UpdateControls(a=rv)
        self.windse_solver.J = 0 # # ALways zero out!!!!
        self.windse_problem.ComputeFunctional()
        self.windse_solver.Solve()
        return self.windse_solver.J

    def eval_QoIGradient(self, mu, xi):
        rv = mu + xi
        self.windse_problem.farm.UpdateControls(a=rv)
        self.windse_solver.J = 0 # ALways zero out!!!!
        self.windse_problem.ComputeFunctional()
        self.windse_solver.Solve()
        self.windse_opt.RecomputeReducedFunctional()
        val = self.windse_opt.Gradient()
        return val

    def plot_eigenvector_on_farm(self, eigenvec, show=True):
        """
        Plotting script for ploting eigenvectors on the windfarm. This is a
        temporary home for this function as it doesnt really belong in here.
        """
        # Get the wind farm bounds
        x = self.windse_problem.farm.x
        y = self.windse_problem.farm.y

        plt.figure()
        # if hasattr(self.dom,"boundary_line"):
        #     plt.plot(*self.dom.boundary_line,c="k")
        # plt.plot(ex_list_x,ex_list_y,c="r")
        p=plt.scatter(x, y, c=eigenvec,cmap="coolwarm",edgecolors=(0, 0, 0, 1))
        # p=plt.scatter(self.x,self.y,c="k",s=70)
        plt.xlim(self.windse_problem.farm.dom.x_range[0],self.windse_problem.farm.dom.x_range[1])
        plt.ylim(self.windse_problem.farm.dom.y_range[0],self.windse_problem.farm.dom.y_range[1])
        clb = plt.colorbar(p)
        clb.ax.set_ylabel('eigenvector value')

        plt.title("Eigenvector Overlay")
        plt.savefig("eigenvector_overlay.pdf", transparent=True)
        if show:
            plt.show()

    def create_solver_object(self):
        if hasattr(self, 'windse_solver'):
            pass
        else:
            self.windse_solver = solve_problem(self.windse_params, self.windse_problem)
            self.windse_opt = windse.Optimizer(self.windse_solver) # You need this to compute gradients


if __name__ == '__main__':

    import copy
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    params_file = "./3x3/2d_wind_farm_induction_factor_opt.yaml"
    n_turbines = 9
    windfarm = WindFarm(n_turbines, params_file, initial_solve=True, mpi_comm=comm)
    og_J = copy.deepcopy(windfarm.windse_solver.J)

    windfarm.plot_eigenvector_on_farm(np.ones(n_turbines))

    rv_val = 0.25*np.ones(n_turbines) # + 0.1*np.random.randn(n_turbines)
    fval = windfarm.eval_QoI(rv_val, np.zeros(n_turbines))
    err = og_J - fval

    grad_val = windfarm.eval_QoIGradient(rv_val, np.zeros(n_turbines))

    # Print it
    print('fval = ', fval)
    print("err = ", err)
    print('grad_val = \n', grad_val)
