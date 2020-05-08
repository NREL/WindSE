.. _params:

The Parameter File
==================

This is a comprehensive list of all the available parameters. The default values are stored within ``default_parameters.yaml`` located in the windse directory of the python source files 


.. contents:: :local:


General Options
---------------

This section is for options about the run itself. The basic format is:: 

    general: 
        name:               <str>
        preappend_datetime: <bool>
        output:             <str list>
        output_folder:      <str> 
        output_type:        <str>
        dolfin_adjoint:     <bool>

+------------------------+-----------------------------------------------------+----------+-----------------------+
| Option                 | Description                                         | Required | Default               |
+========================+=====================================================+==========+=======================+
| ``name``               | Name of the run and the folder in ``output/``.      | no       | "Test"                |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``preappend_datetime`` | Append the date to the output folder name.          | no       | False                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``output``             | | Determines which functions to save.               | no       | ["solution"]          |
|                        | | Select any combination of the following:          |          |                       |
|                        | |   "mesh", "initial_guess", "height",              |          |                       |
|                        | |   "turbine_force", "solution", "debug"            |          |                       |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``output_folder``      | The folder location for all the output              | no       | "output/"             |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``output_type``        | Output format: "pvd" or "xdmf".                     | no       | "pvd"                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+
| ``dolfin_adjoint``     | Required if performing any optimization.            | no       | False                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+



Domain Options
--------------

This section will define all the parameters for the domain::

    domain: 
        type:             <str>
        path:             <str>
        mesh_path:        <str>
        terrain_path:     <str>
        bound_path:       <str>
        filetype:         <str>
        scaled:           <bool>
        ground_reference: <float>
        x_range:          <float list>
        y_range:          <float list>
        z_range:          <float list>
        nx:               <int>
        ny:               <int>
        nz:               <int>
        mesh_type:        <str>
        center:           <float list>
        radius:           <float>
        nt:               <int>
        res:              <int>
        interpolated:     <bool>
        analytic:         <bool>
        gaussian: 
            center:   <float list>
            theta:    <float>
            amp:      <float>
            sigma_x:  <float>
            sigma_y:  <float>
        plane:
            intercept: <float list>
            mx:        <float>
            my:        <float>

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
| ``terrain_path``       | | Location of specific terrain file           | | no               |             | \-          |
|                        | | Default file name: "terrain.txt"            | | "imported"       | ``path``    |             |
|                        | | Note: Only file required by "interpolated"  |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``bound_path``         | | Location of specific boundary marker data   | | no               |             | \-          |
|                        | | Default file name: "boundaries"             | | "imported"       | ``path``    |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``filetype``           | file type for imported mesh: "xml.gz", "h5"   | | no               | "xml.gz"    | \-          |
|                        |                                               | | "imported"       |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``scaled``             | | Scales the domain to km instead of m.       | no                 | False       | \-          |
|                        | | WARNING: extremely experimental!            |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``ground_reference``   | | The height (z coordinate) that is           | no                 | 0.0         | m           |
|                        | | considered ground                           |                    |             |             |
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
| ``interpolated``       | | Indicate if the topography is interpoalted  | | no               |             | \-          |
|                        | | from file or function.                      | | "box"            | False       |             |
|                        |                                               | | "cylinder"       |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``analytic``           | | Indicates if the interpolated function is   | no                 | False       | \-          |
|                        | | analytic or from file.                      |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+

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

+------------------------+-----------------------------------------------+--------------------+---------------+-----------+
| ``plane``              | | If analytic is true, the ground will be     | | "interpolated"   | None          | \-        |
|                        | | represented as a plane                      | | "analytic"       |               |           |
|                        | | Note: requires interpolated and analytic.   |                    |               |           |
+------------------------+-----------------------------------------------+--------------------+---------------+-----------+
| ``intercept``          | The equation of a plane intercept             | no                 | [0.0,0.0,0.0] | m         |
+------------------------+-----------------------------------------------+--------------------+---------------+-----------+
| ``mx``                 | The slope in the x direction                  | yes                | None          | m         |
+------------------------+-----------------------------------------------+--------------------+---------------+-----------+
| ``my``                 | The slope in the y direction                  | yes                | None          | m         |
+------------------------+-----------------------------------------------+--------------------+---------------+-----------+

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
        type:           <str>
        path:           <str>
        display:        <str>
        ex_x:           <float list>
        ex_y:           <float list>
        grid_rows:      <int>
        grid_cols:      <int>
        jitter:         <float>
        numturbs:       <int>
        seed:           <int>
        HH:             <float>
        RD:             <float>
        thickness:      <float>
        yaw:            <float>
        axial:          <float>
        force:          <str>
        turbine_method: <str>
        rpm:            <float>
        read_turb_data: <str>

+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| Option                 | Description                                   | Required (for)     | Default  | Units       |
|                        |                                               |                    |          |             |
+========================+===============================================+====================+==========+=============+
| ``type``               | | Sets the type of farm. Choices:             | yes                | None     | \-          |
|                        | |   "grid", "random", "imported"              |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``path``               | Location of the wind farm text file           | "imported"         | None     | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``display``            | | Displays a plot of the wind farm            | no                 | False    | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``ex_x``               | | The x extents of the farm where turbines    | | "grid"           | None     | m           |
|                        | | can be placed                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``ex_y``               | | The y extents of the farm where turbines    | | "grid"           | None     | m           |
|                        | | can be placed                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``grid_rows``          | The number of turbines in the x direction     | "grid"             | None     | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``grid_cols``          | The number of turbines in the y direction     | "grid"             | None     | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``jitter``             | | Displaces turbines in a random direction    | | no               | 0.0      | m           |
|                        | | by this amount                              | | "grid"           |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``numturbs``           | The total number of turbines                  | "random"           | None     | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``seed``               | | The random seed used to generate/jitter the | | no               | None     | \-          |
|                        | | farm. Useful for repeating random runs      | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``HH``                 | The hub height of the turbine from ground     | | "grid"           | None     | m           |
|                        |                                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``RD``                 | The rotor diameter                            | | "grid"           | None     | m           |
|                        |                                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``thickness``          | The effective thickness of the rotor disk     | | "grid"           | None     | m           |
|                        |                                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``yaw``                | | Determins the yaw of all turbines. Yaw is   | | "grid"           | None     | rad         |
|                        | | relative to the wind inflow direction       | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``axial``              | The axial induction factor                    | | "grid"           | None     | \-          |
|                        |                                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``force``              | | the radial distribution of force            | no                 | "sine"   | \-          |
|                        | | Choices: "sine", "constant"                 |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``turbine_method``     | | determines how the turbine force is built   | no                 | "dolfin" | \-          |
|                        | | Choices: "numpy", "dolfin" , "alm"          |                    |          |             |
|                        | | "numpy"  - builds entirely using arrays,    |                    |          |             |
|                        | |            works best for small farms       |                    |          |             |
|                        | | "dolfin" - uses the FEniCS backend,         |                    |          |             |
|                        | |            robust but potentially slow      |                    |          |             |
|                        | | "alm" - an actuator line method using       |                    |          |             |
|                        | |         numpy array, currently only         |                    |          |             |
|                        | |         support single turbine farms        |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``rpm``                | | sets the revolutions per minute if using    | "alm"              | 10.0     | rev/min     | 
|                        | | the alm turbine method                      |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``read_turb_data``     | | Path to .csv file with chord, lift, and     | no                 | None     | \-          |
|                        | | drag coefficients                           |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+

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
        turbine_type:   <str>
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
|                        | | Choices:                                    |
|                        | | "full" - refines the full mesh              |
|                        | | "box" - refines in a box near the farm      |
|                        | | "cylinder" - cylinder centered at the farm  |
|                        | | "stream" - stream-wise cylinder around farm |
|                        | |            (use for 1 row farms)            |
+------------------------+-----------------------------------------------+
| ``farm_factor``        | | A scaling factor to make the refinement     |
|                        | | area larger or smaller                      |
+------------------------+-----------------------------------------------+
| ``turbine_num``        | Number of turbine refinements                 |
+------------------------+-----------------------------------------------+
| ``turbine_type``       | | The shape of the refinement around turbines |
|                        | | Choices:                                    |
|                        | | "simple" - cylinder around turbine          |
|                        | | "tear" - tear drop shape around turbine     |
|                        | | "wake" - cylinder to capture wake           |
+------------------------+-----------------------------------------------+
| ``turbine_factor``     | | A scaling factor to make the refinement     |
|                        | | area larger or smaller                      |
+------------------------+-----------------------------------------------+
| ``refine_custom``      | | This is a way to define multiple refinements|
|                        | | in a specific order allowing for more       |
|                        | | complex refinement options. Example below   |
+------------------------+-----------------------------------------------+

To use the "refine_custom" option, define a list of lists where each element defines
refinement based on a list of parameters. Example::

    refine_custom: [
        [ "full",     [ ]                                 ],
        [ "full",     [ ]                                 ],
        [ "box",      [ [[-500,500],[-500,500],[0,150]] ] ],
        [ "cylinder", [ [0,0,0], 750, 150 ]               ],
        [ "simple",   [ 100 ]                             ],
        [ "tear",     [ 50, 0.7853 ]                      ]
    ]

For each refinement, the first option indicates how many time this specific
refinement will happen. The second option indicates the type of refinement:
"full", "square", "circle", "farm_circle", "custom". The last option 
indicates the extent of the refinement. 

The example up above will result in five refinements:

    1. Two full refinements
    2. One box refinement bounded by: [[-500,500],[-500,500],[0,150]]
    3. One cylinder centered at origin with radius 750 m and a height of 150 m
    4. One simple turbine refinement with radius 100 m 
    5. One teardrop shaped turbine refinement radius 500 m and rotated by 0.7853 rad

The syntax for each refinement type is::

        [ "full",     [ ]                                                             ]
        [ "box",      [ [[x_min,x_max],[y_min,y_max],[z_min,z_max]], expand_factor ]  ]
        [ "cylinder", [ [c_x,c_y,c_z], radius, height, expand_factor ]                ]
        [ "stream",   [ [c_x,c_y,c_z], radius, length, theta, offset, expand_factor ] ]
        [ "simple",   [ radius, expand_factor ]                                       ]
        [ "tear",     [ radius, theta, expand_factor ]                                ]
        [ "wake",     [ radius, length, theta, expand_factor ]                        ]

Notes::

    * For cylinder, the center is the base of the cylinder
    * For stream, the center is the start of the vertical base and offset indicates the rotation offset
    * For stream, wake, length is the distance center to the downstream end of the cylinder
    * For stream, tear, wake, theta rotates the shape around the center

Function Space Options
----------------------

This section list the function space options::

    function_space:
        type: <str>
        quadrature_degree: <int>
        turbine_space:     <str>
        turbine_degree:    <int>

+------------------------+----------------------------------------------------------+--------------+------------+
| Option                 | Description                                              | Required     | Default    |
|                        |                                                          |              |            |
+========================+==========================================================+==============+============+
| ``type``               | | Sets the type of farm. Choices:                        | yes          | None       |
|                        | |   "linear": P1 elements for both velocity and pressure |              |            |
|                        | |   "taylor_hood": P2 for velocity, P1 for pressure      |              |            |
+------------------------+----------------------------------------------------------+--------------+------------+
| ``quadrature_degree``  | | Sets the quadrature degree for all integration and     | no           | 6          |
|                        | | interpolation for the whole simulation                 |              |            |
+------------------------+----------------------------------------------------------+--------------+------------+
| ``turbine_space``      | | Sets the function space for the turbine. Only needed   | no           | Quadrature |
|                        | | if using "numpy" for ``turbine_method``                |              |            |
|                        | | Choices: "Quadrature", "CG"                            |              |            |
+------------------------+----------------------------------------------------------+--------------+------------+
| ``turbine_degree``     | | The quadrature degree for specifically the turbine     | no           | 6          |
|                        | | force representation. Only works "numpy" method        |              |            |
|                        | | Note: if using Quadrature space, this value must equal |              |            |
|                        | | the ``quadrature_degree``                              |              |            |
+------------------------+----------------------------------------------------------+--------------+------------+



Boundary Condition Options
--------------------------

This section describes the boundary condition options. There are three types
of boundary condtions: inflow, no slip, no stress. By default, inflow is 
prescribed on boundary facing into the wind, no slip on the ground and 
no stress on all other faces. These options describe the inflow boundary
velocity profile.::

    boundary_conditions:
        vel_profile:    <str>
        HH_vel:         <float>
        power:          <float>
        k:              <float>
        turbsim_path    <str>
        inflow_angle:   <float, list>
        boundary_names: <dict>
        boundary_types: <dict>

+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| Option                 | Description                                                                                   | Required     | Default    |
|                        |                                                                                               |              |            |
+========================+===============================================================================================+==============+============+
| ``vel_profile``        | | Sets the velocity profile. Choices:                                                         | yes          | None       |
|                        | |   "uniform": constant velocity of :math:`u_{HH}`                                            |              |            |
|                        | |   "power": a power profile                                                                  |              |            |
|                        | |   "log": log layer profile                                                                  |              |            |
|                        | |   "turbsim": use a turbsim simulation as inflow                                             |              |            |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``HH_vel``             | The velocity at hub height, :math:`u_{HH}`, in m/s.                                           | no           | 8.0        |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``power``              | The power used in the power flow law                                                          | no           | 0.25       |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``k``                  | The constant used in the log layer flow                                                       | no           | 0.4        |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``inflow_angle``       | | Sets the initial inflow angle for the boundary condition. A multiangle solve can be         | no           | None       |
|                        | | indicated by setting this value to a list with values: [start, stop, n] where the solver    |              |            |
|                        | | will perform n solves, sweeping uniformly through the start and stop angles. The number of  |              |            |
|                        | | solves, n, can also be defined in the solver parameters.                                    |              |            |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``turbsim_path``       | The location of turbsim profiles used as inflow boundary conditions                           | | yes        | None       |
|                        |                                                                                               | | "turbsim"  |            |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``boundary_names``     | A dictionary used to identify the boundaries                                                  | no           | See Below  |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+
| ``boundary_types``     | A dictionary for defining boundary conditions                                                 | no           | See Below  |
+------------------------+-----------------------------------------------------------------------------------------------+--------------+------------+

..
    of :math:`u_x=u_{max} \left( \frac{z-z_0}{z_1-z_0} \right)^{p}`

If you are importing a mesh or want more control over boundary conditions, you can specify the boundary markers using ``names`` and ``types``.
The default for these two are

Rectangular Mesh::

    boundary_condition:
        boundary_names: 
            east:  1
            north: 2
            west:  3
            south: 4
        boundary_types: 
            inflow:    ["west","north","south"]
            no_stress: ["east"]

Box Mesh::

    boundary_condition:
        boundary_names: 
            east:   1
            north:  2
            west:   3
            south:  4
            bottom: 5
            top:    6
        boundary_types: 
            inflow:    ["west","north","south"]
            free_slip: ["top"]
            no_slip:   ["bottom"]
            no_stress: ["east"]

Circle Mesh::

    boundary_condition:
        boundary_names: 
            outflow: 7
            inflow:  8
        boundary_types: 
            inflow:    ["inflow"]
            no_stress: ["outflow"]

Cylinder Mesh::

    boundary_condition:
        boundary_names: 
            outflow: 5
            inflow:  6
            bottom:  7
            top:     8
        boundary_types: 
            inflow:    ["inflow"]
            free_slip: ["top"]
            no_slip:   ["bottom"]
            no_stress: ["outflow"]

These defaults correspond to an inflow wind direction from West to East.

When marking a rectangular/box domains, from a top-down perspective, start from 
the boundary in the positive x direction and go counter clockwise, the boundary 
names are: "easy", "north", "west", "south". Additionally, in 3D there are also
"top" and "bottom". For a circular/cylinder domains, the boundary names are
"inflow" and "outflow". Likewise, in 3D there are also "top" and "bottom". 
Additionally, you can change the ``boundary_types`` if using one of the built 
in domain types. This way you can customize the boundary conditions without 
importing a whole new mesh.

Problem Options
---------------

This section describes the problem options::

    problem:
        type:          <str>
        use_25d_model: <bool>
        viscosity:     <float>
        lmax:          <float>

+------------------------+--------------------------------------------------------------+--------------+---------+
| Option                 | Description                                                  | Required     | Default |
|                        |                                                              |              |         |
+========================+==============================================================+==============+=========+
| ``type``               | | Sets the variational form use. Choices:                    | yes          | None    |
|                        | |   "taylor_hood": Standard RANS formulation                 |              |         |
|                        | |   "stabilized": Adds a term to stabilize P1xP1 formulations|              |         |
+------------------------+--------------------------------------------------------------+--------------+---------+
| ``viscosity``          | Kinematic Viscosity                                          | no           | 0.1     |
|                        |                                                              |              |         |
+------------------------+--------------------------------------------------------------+--------------+---------+
| ``lmax``               | Turbulence length scale                                      | no           | 15.0    |
|                        |                                                              |              |         |
+------------------------+--------------------------------------------------------------+--------------+---------+
| ``use_25d_model``      | | Option to enable a small amount of compressibility to mimic| no           | False   |
|                        | | the effect of a 3D, out-of-plane flow solution in a 2D     | "2D only"    |         |
|                        | | model.                                                     |              |         |
+------------------------+--------------------------------------------------------------+--------------+---------+




Solver Options
--------------

This section lists the solver options::

    solver:
        type:              <str>
        final_time:        <float>
        record_time:       <str, float>
        save_interval:     <float>
        num_wind_angles:   <int>
        endpoint:          <bool>
        velocity_path:     <str>
        save_power:        <bool>
        nonlinear_solver:  <str>
        newton_relaxation: <float>

+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| Option                 | Description                                                    | Required (for)      | Default             |
|                        |                                                                |                     |                     |
+========================+================================================================+=====================+=====================+
| ``type``               | | Sets the solver type. Choices:                               | yes                 | None                |
|                        | |   "steady": solves for the steady state solution             |                     |                     |
|                        | |   "unsteady": solves for a time varying solution             |                     |                     |
|                        | |   "multiangle": iterates through inflow angles               |                     |                     |
|                        | |                 uses ``inflow_angle`` or [0, :math:`2\pi`]   |                     |                     |
|                        | |   "imported_inflow": runs multiple steady solves with        |                     |                     |
|                        | |                      imported list of inflow conditions      |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``final_time``         | The final time for an unsteady simulation                      | | no                | 1.0 s               |
|                        |                                                                | | "unsteady"        |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``record_time``        | | The time when we start recording the objective function      | | no                | None                |
|                        | | using a weighted average.                                    | | "unsteady"        |                     |
|                        | | Choices:                                                     |                     |                     |
|                        | |          "computed": Sets the time based on inflow speed     |                     |                     |
|                        | |          "last": Only computes at the final time step        |                     |                     |
|                        | |          <float>: use this specific time to start recording  |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``save_interval``      | The amount of time between saving output fields                | | no                | 1.0 s               |
|                        |                                                                | | "unsteady"        |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``num_wind_angles``    | Sets the number of angles. can also be set in ``inflow_angle`` | | no                | 1                   |
|                        |                                                                | | "multiangle"      |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``endpoint``           | Should the final inflow angle be simulated                     | | no                | False               |
|                        |                                                                | | "multiangle"      |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``velocity_path``      | The location of a list of inflow conditions                    | | yes               |                     |
|                        |                                                                | | "imported_inflow" |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``save_power``         | Save the power for each turbine to a text file in              | no                  | False               |
|                        | output/``name``/data/                                          |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``nonlinear_solver``   | | Specify the nonlinear solver type. Choices:                  | no                  | "snes"              |
|                        | |   "newton": uses the standard newton solver                  |                     |                     |
|                        | |   "snes": PETSc SNES solver                                  |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``newton_relaxation``  | Set the relaxation parameter if using newton solver            | | no                | 1.0                 |
|                        |                                                                | | "newton"          |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+

The "multiangle" solver uses the steady solver to solve the RANS formulation.
Currently, the "multiangle" solver does not support imported domains. 


Optimization Options
--------------------

This section lists the optimization options. If you are planning on doing
optimization make sure to set ``dolfin_adjoint`` to True.::

    optimization:
        control_types:  <str list>
        layout_bounds:  <float list>
        objective_type: <str>
        wake_RD:        <int>
        min_total:      <int>
        taylor_test:    <bool>
        optimize:       <bool>
        gradient:       <bool>

+------------------------+----------------------------------------------------------+-----------------+--------------+
| Option                 | Description                                              | Required        | Default      |
|                        |                                                          |                 |              |
+========================+==========================================================+=================+==============+
| ``control_types``      | | Sets the parameters to optimize. Choose Any:           | yes             | None         |
|                        | |   "yaw", "axial", "layout"                             |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``layout_bounds``      | The bounding box for the layout optimization             | no              | wind_farm    |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``objective_type``     | | Sets the objective function for optimization           | no              | power        |
|                        | | Choices:                                               |                 |              |
|                        | |   "power": simple power calculation                    |                 |              |
|                        | |   "2d_power": power calculation optimized for 2D runs  |                 |              |
|                        | |   "wake_deflection": metric for measuring wake movement|                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``wake_RD``            | | number of rotor diameters downstream where the wake is | no              | 5            |
|                        | | measured                                               | wake_deflection |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``min_total``          | | number of times the average wake deflection reaches a  | no              | 0            |
|                        | | before the unsteady simulation is stopped. use 0 to    | wake_deflection |              |
|                        | | run the full simulation                                |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``taylor_test``        | | Performs a test to check the derivatives. Good         | no              | False        |
|                        | | results have a convergence rate around 2.0             |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``optimize``           | | Optimize the given controls using the power output as  | no              | False        |
|                        | | the objective function using SLSQP from scipy.         |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``gradient``           | | returns the gradient values of the objective with      | no              | False        |
|                        | | respect to the controls                                |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
