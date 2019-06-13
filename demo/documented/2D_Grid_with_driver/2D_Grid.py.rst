
.. _demo_2d_grid:

Gridded Wind Farm on a Rectangular Domain
=========================================

This demonstration will show how to set up a 2D rectangular mesh with a
wind farm consisting of a 36 turbines laid out in a 6x6 grid. Running this 
demo requires two files:
    
    * Parameter File: :download:`params.yaml`
    * Driver File: :download:`2D_Grid.py`


Setting up the parameters file:
-------------------------------

To use WindSE, each simulation needs a driver file and a parameters file.
The parameter file is the main way to customize a simulation. The driver 
file uses the options specified in the parameters file to run the simulation.
Ideally, multiple simulations can use a single driver file and multiple 
parameter files. 

The parameter file is formated as a `yaml <https://yaml.org/>` structure and
requires `pyyaml <https://pyyaml.org/>` to be read. The driver file is 
written in python. 

The parameter file is broken up into several sections: general, domain, 
boundaries, and wind_farm. For this demo we will will not need to use the 
boundaries section. 

The full parameter file can be found here: :download:`params.yaml`

.. include:: parameters.rst

Creating the driver code:
-------------------------

The full driver file can be found here: :download:`2D_Grid.py` First, 
we start off with the import statements::

    import windse

Since, this driver is not that fancy, we only need to import :py:mod:`windse`
Next we need to initialize WindSE with the parameters::

    windse.initialize("params.yaml")

which are located in the root directory. However, we could have placed the 
parameters files any where with any name as long as we direct WindSE to that
particular file.

Now we can generate the domain and wind farm::

    dom = windse.RectangleDomain()
    farm = windse.GridWindFarm(dom)

We can inspect the wind farm by running::

    farm.Plot(True)

This results in a wind farm that looks like this:

.. figure:: wind_farm.png
   :scale: 75 %

Alternatively, we could have use ``False`` to generate and save the plot, 
but not display it. This is useful for running batch test or on a HPC. We 
could also manually save the mesh using ``dom.Save()``, but since we 
specified the mesh as an output in the parameters file, this will be done
automatically when we solve.  

Next we need to build the finite element function space::

    fs = windse.TaylorHoodFunctionSpace2D(dom)

For this problem we are going to use Taylor-Hood elements, which are
comprised of 2nd order Lagrange elements for velocity and 1st order elements
for pressure. 

The last step before assembling the problem is to define the boundary
conditions::

    bc = windse.UniformInflow2D(dom,fs)

This problem has uniform inflow from the left and includes the top, left and
bottom boundaries. The right boundary is our outflow and has a no-stress
boundary condition.

Finally, it's time to assemble everything and solve::
    
    problem = windse.TaylorHoodProblem2D(dom,farm,fs,bc)
    solver = windse.SteadySolver(problem)
    solver.Solve()

Running ``solver.Solve()`` will save all the inputs according to the 
parameters file, solve the problem, and save the solution. If everything 
went smoothly, the solution for wind speed should be:


.. image:: solution.png
   :scale: 75 %