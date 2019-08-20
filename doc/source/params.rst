.. _params:

The Parameter File
==================

.. contents:: :local:


General Options
---------------

This section is for options about the run itself. The basic format is:: 

    general: 
        name:               <str>
        preappend_datetime: <bool>
        output:             <str list>
        output_type:        <str>
        dolfin_adjoint:     <bool>

+------------------------+-----------------------------------------------------+----------+-----------------------+
| Option                 | Description                                         | Required | Default               |
+========================+=====================================================+==========+=======================+
| ``name``               | Name of the run and the folder in ``output/``.      | no       | Test                  |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``preappend_datetime`` | Append the date to the output folder name.          | no       | False                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``output``             | | Determines which functions to save.               | no       | ["solution"]          |
|                        | | Select any combination of the following:          |          |                       |
|                        | |   "mesh", "initial_guess", "height",              |          |                       |
|                        | |   "turbine_force", "solution", "debug"            |          |                       |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``output_type``        | Output format: "pvd" or "xdmf".                     | no       | "pvd"                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``dolfin_adjoint``     | Required if performing any optimization.            | no       | False                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+



Domain Options
--------------

This section will define all the parameters for the domain::

    domain: 
        type:         <str>
        path:         <str>
        mesh_path:    <str>
        typo_path:    <str>
        bound_path:   <str>
        filetype:     <str>
        x_range:      <float list>
        y_range:      <float list>
        z_range:      <float list>
        nx:           <int>
        ny:           <int>
        nz:           <int>
        mesh_type:    <str>
        center:       <float list>
        radius:       <float>
        nt:           <int>
        res:          <int>
        interpolated: <bool>
        analytic:     <bool>
        gaussian: 
            center:   <float list>
            amp:      <float>
            sigma_x:  <float>
            sigma_y:  <float>
            theta:    <float>

+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| Option                 | Description                                   | Required (for)     | Default     | Units       |
+========================+===============================================+====================+=============+=============+
| ``type``               | | Sets the shape/dimension of the mesh.       | yes                | None        | \-          |
|                        | | Choices:                                    |                    |             |             |
|                        | |   "rectangle", "box", "cylinder", "circle"  |                    |             |             |
|                        | |   "imported", "interpolated"                |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``path``               | Folder of the mesh data to import             | | yes or ``*_path``|             | \-          |
|                        |                                               | | "imported"       | ``*_path``  |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``mesh_path``          | | Location of specific mesh file              | | no               |             | \-          |
|                        | | Default file name: "mesh"                   | | "imported"       | ``path``    |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``typo_path``          | | Location of specific topology file          | | no               |             | \-          |
|                        | | Default file name: "topology.txt"           | | "imported"       | ``path``    |             |
|                        | | Note: Only file required by "interpolated"  |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``bound_path``         | | Location of specific boundary marker data   | | no               |             | \-          |
|                        | | Default file name: "boundaries"             | | "imported"       | ``path``    |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``filetype``           | file type for imported mesh: "xml.gz", "h5"   | | no               | "xml.gz"    | \-          |
|                        |                                               | | "imported"       |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``x_range``            | List of two floats defining the x range       | | "rectangle"      | None        | m           |
|                        |                                               | | "box"            |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``y_range``            | List of two floats defining the y range       | | "rectangle"      | None        | m           |
|                        |                                               | | "box"            |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``z_range``            | List of two floats defining the z range       | | "box"            | None        | m           |
|                        |                                               | | "cylinder"       |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``nx``                 | The number of nodes in the x direction        | | "rectangle"      | None        | \-          |
|                        |                                               | | "box"            |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``ny``                 | The number of nodes in the x direction        | | "rectangle"      | None        | \-          |
|                        |                                               | | "box"            |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``nz``                 | The number of nodes in the x direction        | | "box"            | None        | \-          |
|                        |                                               | | "cylinder"       |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``mesh_type``          | | The meshing type when generating a          | | "cylinder"       | "mshr"      | \-          |
|                        | | cylindric domain.                           | | "circle"         |             |             |
|                        | | Choices:                                    |                    |             |             |
|                        | |   "mshr", "elliptic", "squircular",         |                    |             |             |
|                        | |   "stretch"                                 |                    |             |             |
|                        | | Note: ``nz`` doesn't work with "mshr"       |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``center``             | A 2D list indicating the center of the base   | | "cylinder"       | None        | m           |
|                        |                                               | | "circle"         |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``radius``             | The radius of the cylinder                    | | "cylinder"       | None        | m           |
|                        |                                               | | "circle"         |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``nt``                 | | The number of radial segments to            | | "cylinder"       | None        | \-          |
|                        | | approximate the cylinder                    | | "circle"         |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``res``                | | The resolution of the mesh. It should be    | | "cylinder"       | None        | \-          |
|                        | | less than ``nt``.                           | | "circle"         |             |             |
|                        | | Note: ``res`` only works with "mshr"        |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``interpolated``       | |Indicate if the topography is interpoalted   | | no               |             | \-          |
|                        | |from file or function.                       | | "box"            | False       |             |
|                        |                                               | | "cylinder"       |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``analytic``           | |Indicates if the interpolated function is    | "interpolated"     |             | \-          |
|                        | |analytic.                                    |                    | False       |             |
|                        |                                               |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``gaussian``           | | If analytic is true, a Gaussian hill will   | | "interpolated"   | None        | \-          |
|                        | | be created using the following parameters.  | | "analytic"       |             |             |
|                        | | Note: requires interpolated and analytic.   |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``center``             | The center point of the gaussian hill.        | no                 | [0.0,0.0]   | m           |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``amp``                | The amplitude of the hill.                    | yes                | None        | m           |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``sigma_x``            | The extent of the hill in the x direction.    | yes                | None        | m           |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``sigma_y``            | The extent of the hill in the y direction.    | yes                | None        | m           |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``theta``              | The rotation of the hill.                     | no                 | 0.0         | rad         |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+

To import a domain, three files are required: 

* mesh.xml.gz - this contains the mesh in a format dolfin can handle
* boundaries.xml.gz - this contains the facet markers that define where the boundaries are
* topology.txt - this contains the data for the ground topology. 

The topology file assumes that the coordinates are from a uniform mesh.
It contains three column: x, y, z. The x and y columns contain 
just the unique values. The z column contains the ground values
for every combination of x and y. The first row must be the number
of points in the x and y direction. Here is an example for z=x+y/10::

            3 3 9
            0 0 0.0
            1 1 0.1
            2 2 0.2
                1.0
                1.1
                1.2
                2.0
                2.1
                2.2

Note: If using "h5" file format, the mesh and boundary will be in one file.



Wind Farm Options
-----------------

This section will define all the parameters for the wind farm::

    wind_farm: 
        type:      <str>
        path:      <str>
        display:   <str>
        ex_x:      <float list>
        ex_y:      <float list>
        grid_rows: <int>
        grid_cols: <int>
        jitter:    <float>
        numturbs:  <int>
        seed:      <int>
        HH:        <float>
        RD:        <float>
        thickness: <float>
        yaw:       <float>
        axial:     <float>

+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| Option                 | Description                                   | Required (for)     | Default | Units       |
|                        |                                               |                    |         |             |
+========================+===============================================+====================+=========+=============+
| ``type``               | | Sets the type of farm. Choices:             | yes                | None    | \-          |
|                        | |   "grid", "random", "imported"              |                    |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``path``               | Location of the wind farm text file           | "imported"         | None    | \-          |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``display``            | | Displays a plot of the wind farm            | no                 | False   | \-          |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``ex_x``               | | The x extents of the farm where turbines    | | "grid"           | None    | m           |
|                        | | can be placed                               | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``ex_y``               | | The y extents of the farm where turbines    | | "grid"           | None    | m           |
|                        | | can be placed                               | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``grid_rows``          | The number of turbines in the x direction     | "grid"             | None    | \-          |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``grid_cols``          | The number of turbines in the y direction     | "grid"             | None    | \-          |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``jitter``             | | Displaces turbines in a random direction    | | no               | 0.0     | m           |
|                        | | by this amount                              | | "grid"           |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``numturbs``           | The total number of turbines                  | "random"           | None    | \-          |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``seed``               | | The random seed used to generate the farm.  | | no               | None    | \-          |
|                        | | Useful for repeating random runs            | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``HH``                 | The hub height of the turbine from ground     | | "grid"           | None    | m           |
|                        |                                               | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``RD``                 | The rotor diameter                            | | "grid"           | None    | m           |
|                        |                                               | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``thickness``          | The effective thickness of the rotor disk     | | "grid"           | None    | m           |
|                        |                                               | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``yaw``                | | Determins the yaw of all turbines. Yaw is   | | "grid"           | None    | rad         |
|                        | | relative to the wind inflow direction       | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+
| ``axial``              | The axial induction factor                    | | "grid"           | None    | \-          |
|                        |                                               | | "random"         |         |             |
+------------------------+-----------------------------------------------+--------------------+---------+-------------+

To import a wind farm, create a .txt file with this formatting::

    #    x      y     HH    Yaw   Diameter Thickness Axial_Induction
    200.00 0.0000 80.000  0.000      126.0      10.5            0.33
    800.00 0.0000 80.000  0.000      126.0      10.5            0.33

The first row isn't necessary. Each row defines a different turbine.



Refinement Options
------------------

This section describes the options for refinement
The domain created with the previous options can be refined in special
ways to maximize the efficiency of the number DOFs. None of these options
are required. There are three types of mesh manipulation: warp, farm refine,
turbine refine. Warp shifts more cell towards the ground, refining the farm
refines within the farm extents, and refining the turbines refines within
the rotor diameter of a turbine. When choosing to warp, a "smooth" warp will 
shift the cells smoothly towards the ground based on the strength. A "split"
warp will attempt to create two regions, a high density region near the 
ground and a low density region near the top

The options are::

    refine:
        warp_type:      <str>
        warp_strength:  <float>
        warp_percent:   <float>
        warp_height:    <float>
        farm_num:       <int>
        farm_type:      <str>
        farm_factor:    <float>
        turbine_num:    <int>
        turbine_factor: <float>
        refine_custom:  <list list>

+------------------------+-----------------------------------------------+
| Option                 | Description                                   |
+========================+===============================================+
| ``warp_type``          | | Choose to warp the mesh to place more cells |
|                        | | near the ground. Choices:                   |
|                        | |   "smooth", "split"                         |
+------------------------+-----------------------------------------------+
| ``warp_strength``      | | The higher the strength the more cells      |
|                        | | moved towards the ground. Requires: "smooth"|
+------------------------+-----------------------------------------------+
| ``warp_percent``       | | The percent of the cell moved below the     |
|                        | | warp height. Requires: "split"              |
+------------------------+-----------------------------------------------+
| ``warp_height``        | | The height the cell are moved below         |
|                        | | Requires: "split"                           |
+------------------------+-----------------------------------------------+
| ``farm_num``           | Number of farm refinements                    |
+------------------------+-----------------------------------------------+
| ``farm_type``          | | The shape of the refinement around the farm |
|                        | | Choices: "square", "circle", "farm_circle"  |
|                        | | "square" and "circle" are centered at the   |
|                        | | center of the domain, whereas, "farm_circle"|
|                        | | is centered at the farm                     |
+------------------------+-----------------------------------------------+
| ``farm_factor``        | | A scaling factor to make the refinement     |
|                        | | area larger or smaller                      |
+------------------------+-----------------------------------------------+
| ``turbine_num``        | Number of turbine refinements                 |
+------------------------+-----------------------------------------------+
| ``turbine_factor``     | | A scaling factor to make the refinement     |
|                        | | area larger or smaller                      |
+------------------------+-----------------------------------------------+
| ``refine_custom``      | | This is a way to define multiple refinements|
|                        | | in a specific order allowing for more       |
|                        | | complex refinement options. Example below   |
+------------------------+-----------------------------------------------+

To use the "refine_custom" option create a list of where each element defines
refinement based on a list of parameters. Example::

    refine_custom: [
        [2,full],
        [1, custom, [[-1200,1200],[-1200,1200],[0,140]]],
        [1, circle, 1020]
        [1, farm_circle, 1020]
        [1, square, 500]
    ] 

For each refinement, the first option indicates how many time this specific
refinement will happen. The second option indicates the type of refinement:
"full", "square", "circle", "farm_circle", "custom". The last option 
indicates the extent of the refinement. 

The example up above will result in five refinements:

    1. Two full refinements
    2. One custom refinement bounded by the box: [[-1200,1200],[-1200,1200],[0,140]]
    3. One circle centered at the center of the domain with radius 1020 m
    4. One circle centered at the center of the farm with radius 1020 m 
    5. One square centered at the center of the farm with side length 500 m


Function Space Options
----------------------

This section list the function space options::

    function_space:
        type: <str>

+------------------------+----------------------------------------------------------+--------------+---------+
| Option                 | Description                                              | Required     | Default |
|                        |                                                          |              |         |
+========================+==========================================================+==============+=========+
| ``type``               | | Sets the type of farm. Choices:                        | yes          | None    |
|                        | |   "linear": P1 elements for both velocity and pressure |              |         |
|                        | |   "taylor_hood": P2 for velocity, P1 for pressure      |              |         |
+------------------------+----------------------------------------------------------+--------------+---------+



Boundary Condition Options
--------------------------

This section describes the boundary condition options. There are three types
of boundary condtions: inflow, no slip, no stress. By default, inflow is 
prescribed on boundary facing into the wind, no slip on the ground and 
no stress on all other faces. These options describe the inflow boundary
velocity profile.::

    boundary_condition:
        boundary_names: <dict>
        boundary_types: <dict>
        vel_profile:    <str>
        HH_vel:         <float>
        power:          <float>
        k:              <float>

+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| Option                 | Description                                                                                   | Required     | Default    |
|                        |                                                                                               |              |            |
+========================+===============================================================================================+==============+============+
| ``boundary_names``     | A dictionary used to identify the boundaries                                                  | no           | See Below  |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``boundary_types``     | A dictionary for defining boundary conditions                                                 | no           | See Below  |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``vel_profile``        | | Sets the velocity profile. Choices:                                                         | yes          | None       |
|                        | |   "uniform": constant velocity of :math:`u_{HH}`                                            |              |            |
|                        | |   "power": a power profile                                                                  |              |            |
|                        | |   "log": log layer profile                                                                  |              |            |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``HH_vel``             | The velocity at hub height, :math:`u_{HH}`, in m/s.                                           | no           | 8.0        |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``power``              | The power used in the power flow law                                                          | no           | 0.25       |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``k``                  | The constant used in the log layer flow                                                       | no           | 0.4        |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+

..
    of :math:`u_x=u_{max} \left( \frac{z-z_0}{z_1-z_0} \right)^{p}`

If you are importing a mesh, you can specify the boundary markers using ``names`` and ``types``.
The default for these two are
Rectangular Mesh::

    boundary_condition:
        boundary_names: 
            front: 1
            back: 2
            left: 3
            right: 4
        boundary_types: 
            inflow: ["front","left","right"],
            no_stress: ["back"]

Box Mesh::

    boundary_condition:
        boundary_names: 
            top: 1
            bottom: 2
            front: 3
            back: 4
            left: 5
            right: 6
        boundary_types: 
            inflow: ["top","front","left","right"],
            no_slip:   ["bottom"]
            no_stress: ["back"]

Cylinder/Interpolated Mesh::

    boundary_condition:
        boundary_names: 
            inflow: 1
            outflow: 2
            top: 3
            bottom: 4
        boundary_types: 
            inflow: ["inflow","top"],
            no_slip:   ["bottom"]
            no_stress: ["outflow"]

These defaults corrispond to an inflow wind direction from West to East.
Feel free to mimic these defaults when creating the boundary input file. 
Alternatively, you can name you boundaries whatever you want as long as you
set up the corresponding ``boundary_types``. Additionally, you can set 
change the ``boundary_types`` if using one of the built in domain types. 
This way you can customize the boundary conditions without importing a whole
new mesh.



Problem Options
---------------

This section describes the problem options::

    problem:
        type: <str>

+------------------------+--------------------------------------------------------------+--------------+---------+
| Option                 | Description                                                  | Required     | Default |
|                        |                                                              |              |         |
+========================+==============================================================+==============+=========+
| ``type``               | | Sets the variational form use. Choices:                    | yes          | None    |
|                        | |   "taylor_hood": Standard RANS formulation                 |              |         |
|                        | |   "stabilized": Adds a term to stabilize P1xP1 formulations|              |         |
+------------------------+--------------------------------------------------------------+--------------+---------+




Solver Options
--------------

This section lists the solver options::

    solver:
        type:             <str>
        wind_range:       <float list>
        endpoint:         <bool>
        num_wind_angles:  <int>

+------------------------+----------------------------------------------------------+-------------------+---------------------+
| Option                 | Description                                              | Required (for)    | Default             |
|                        |                                                          |                   |                     |
+========================+==========================================================+===================+=====================+
| ``type``               | | Sets the solver type. Choices:                         | yes               | None                |
|                        | |   "steady": solves for the steady state solution       |                   |                     |
|                        | |   "multiangle": iterates through inflow angles         |                   |                     |
+------------------------+----------------------------------------------------------+-------------------+---------------------+
| ``wind_range``         | The start and end angles to sweep over                   | | no              | [0.0, :math:`2\pi`] |
|                        |                                                          | | "multiangle"    |                     |
+------------------------+----------------------------------------------------------+-------------------+---------------------+
| ``endpoint``           | Should the end point be included in the run              | | no              | False               |
|                        |                                                          | | "multiangle"    |                     |
+------------------------+----------------------------------------------------------+-------------------+---------------------+
| ``num_wind_angles``    | Sets the number of angles                                | | yes             | None                |
|                        |                                                          | | "multiangle"    |                     |
+------------------------+----------------------------------------------------------+-------------------+---------------------+

The "multiangle" solver uses the steady solver to solve the RANS formulation.
Currently, the "multiangle" solver does not support imported domains. 


Optimization Options
--------------------

This section lists the optimization options. If you are planning on doing
optimization make sure to set ``dolfin_adjoint`` to True.::

    optimization:
        controls:     <str list>
        layout_bounds <float list>
        taylor_test:  <bool>
        optimize:     <bool>

+------------------------+----------------------------------------------------------+-------------+--------------+
| Option                 | Description                                              | Required    | Default      |
|                        |                                                          |             |              |
+========================+==========================================================+=============+==============+
| ``controls``           | | Sets the parameters to optimize. Choose Any:           | yes         | None         |
|                        | |   "yaw", "axial", "layout"                             |             |              |
+------------------------+----------------------------------------------------------+-------------+--------------+
| ``taylor_test``        | | Performs a test to check the derivatives. Good         | no          | True         |
|                        | | results have a convergence rate around 2.0             |             |              |
+------------------------+----------------------------------------------------------+-------------+--------------+
| ``optimize``           | | Optimize the given controls using the power output as  | no          | True         |
|                        | | the objective function using SLSQP from scipy.         |             |              |
+------------------------+----------------------------------------------------------+-------------+--------------+
| ``layout_bounds``      | The bounding box for the layout optimization             | no          | wind_farm_ex |
+------------------------+----------------------------------------------------------+-------------+--------------+
