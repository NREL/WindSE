
.. _parameter_file:

Setting up general options:
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The general options are those that will effect the entire run and usually
specify how to handle i/o. for this demo the general parameters are::

  general: 
    name: "2D"
    preappend_datetime: false
    output: ["mesh","initial_guess","turbine_force","solution"]
    output_type: "xdmf"

The ``name`` parameter determines the naming structure for the output 
folders. usually the output folder is ``output/<name>/``. This is the only 
required options.

Setting ``preappend_datetime`` to ``true`` will append the ``name`` with a 
datetime stamp. This is useful when running multiple simulation as they will
be organized by date. The default option for this is ``false``

The ``outputs`` is a list of function that will be saved when 
``solver.Solve()`` is called. These strings can be in any combination:
  
  * ``mesh``: saves the mesh and boundary markers
  * ``initial_guess``: saves the initial velocity and pressure used by the Newton iteration
  * ``height``: saves a function indicating the terrain height and depth
  * ``turbine_force``: saves the function that is used to represent the turbines
  * ``solution``: saves the velocity and pressure after a solve

By default, the only output is ``solution``.

Finally, the ``output_type`` is the file format for the saved function. 
Currently WindSE supports ``xdmf`` and ``pvd`` with the latter being the 
default. However, the mesh files are always saved in the pvd format.

Setting up the domain:
^^^^^^^^^^^^^^^^^^^^^^

Next we need to set the parameters for the domain::

  domain: 
    #                      # Description           | Units
    x_range: [-2500, 2500] # x-range of the domain | m
    y_range: [-2500, 2500] # y-range of the domain | m
    nx: 200                # Number of x-nodes     | -
    ny: 200                # Number of y-nodes     | -

This will create a mesh that has 200 nodes in the x-direction and 200 nodes
in the y-direction. The mesh will be a rectangle with side lengths of 5000 m
and centered at (0,0). 

Setting up the wind farm:
^^^^^^^^^^^^^^^^^^^^^^^^^

The last step for this demo is to set up the wind farm::

  wind_farm: 
    #                  # Description              | Units
    ex_x: [-1800,1800] # x-extent of the farm     | m
    ex_y: [-1800,1800] # y-extent of the farm     | m
    grid_rows: 6       # Number of rows           | -
    grid_cols: 6       # Number of columns        | -
    yaw: 0             # Yaw                      | rads
    axial: 0.33        # Axial Induction          | -
    HH: 90             # Hub Height               | m
    RD: 126            # Turbine Diameter         | m
    thickness: 10      # Effective Thickness      | m

This will produce a 6 by 6 grid evenly spaced in an area of 
[-1800,1800] X [-1800,1800]. Note that ``ex_x`` X ``ex_y`` is the extent of the
farm and should be a subset of the domain ranges. The extent accounts for 
the rotor diameter to ensure all turbines including the rotors are located
within the extents. The rest of the parameters determine the physical 
properties of the turbines:

  * ``yaw``: The yaw of the turbines where 0 is perpendicular to inflow.
  * ``axial``: The axial induction
  * ``HH``: The hub height relative to the ground
  * ``RD``: The rotor diameter
  * ``thickness``: The effective thickness of the rotor used for calculating the turbine force

Other Required Parameters:
^^^^^^^^^^^^^^^^^^^^^^^^^^

Additionally, we need to specify a few parameters that are required for some checks.
These options are not actually used within the custom driver::

  problem:
    type: taylor-hood

  solver:
    type: steady