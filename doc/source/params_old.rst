.. _params:

Legacy Parameter File
=====================

This is a comprehensive list of all the available parameters. The default values are stored within ``default_parameters.yaml`` located in the windse directory of the python source files. 


.. contents:: :local:

Adding a New Parameter
----------------------

To add a new parameter, first add an entry in the ``default_parameters.yaml`` file under one of the major sections with a unique name. The value of that entry will then be an attribute of the class associated with that section. For example: adding the entry ``test_param: 42`` under the ``wind_farm`` section will add the attribute ``self.test_param`` to the ``GenericWindFarm`` class with the default value of 42.


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
        debug_mode:         <bool>

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
| ``debug_mode``         | Saves "tagged_output.yaml" full of debug data       | no       | False                 |
+------------------------+-----------------------------------------------------+----------+-----------------------+



Domain Options
--------------

This section will define all the parameters for the domain::

    domain: 
        type:                <str>
        path:                <str>
        mesh_path:           <str>
        terrain_path:        <str>
        bound_path:          <str>
        filetype:            <str>
        scaled:              <bool>
        ground_reference:    <float>
        streamwise_periodic: <bool>
        spanwise_periodic:   <bool>
        x_range:             <float list>
        y_range:             <float list>
        z_range:             <float list>
        nx:                  <int>
        ny:                  <int>
        nz:                  <int>
        mesh_type:           <str>
        center:              <float list>
        radius:              <float>
        nt:                  <int>
        res:                 <int>
        interpolated:        <bool>
        analytic:            <bool>
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
| ``streamwise_periodic``| | Sets periodic boundary condition in the x   | no                 | False       | \-          |
|                        | | direction (NOT FULLY IMPLEMENTED)           |                    |             |             |
+------------------------+-----------------------------------------------+--------------------+-------------+-------------+
| ``spanwise_periodic``  | | Sets periodic boundary condition in the y   | no                 | False       | \-          |
|                        | | direction (NOT FULLY IMPLEMENTED)           |                    |             |             |
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
| ``interpolated``       | | Indicate if the topography is interpolated  | | no               |             | \-          |
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
        type:               <str>
        path:               <str>
        display:            <str>
        ex_x:               <float list>
        ex_y:               <float list>
        x_spacing:          <float>
        y_spacing:          <float>
        x_shear:            <float>
        y_shear:            <float>
        min_sep_dist:       <float>
        grid_rows:          <int>
        grid_cols:          <int>
        jitter:             <float>
        numturbs:           <int>
        seed:               <int>

+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| Option                 | Description                                   | Required (for)     | Default  | Units       |
|                        |                                               |                    |          |             |
+========================+===============================================+====================+==========+=============+
| ``type``               | | Sets the type of farm. Choices:             | yes                | None     | \-          |
|                        | |   "grid", "random", "imported", "empty"     |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``path``               | Location of the wind farm csv file            | "imported"         | None     | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``display``            | | Displays a plot of the wind farm            | no                 | False    | \-          |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``ex_x``               | | The x extents of the farm where turbines    | | "grid"           | None     | m           |
|                        | | can be placed                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``ex_y``               | | The y extents of the farm where turbines    | | "grid"           | None     | m           |
|                        | | can be placed                               | | "random"         |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``x_spacing``          | | Alternative method for defining grid farm   | "grid"             | None     | m           |
|                        | | x distance between turbines                 |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``y_spacing``          | | Alternative method for defining grid farm   | "grid"             | None     | m           |
|                        | | y distance between turbines                 |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``x_shear``            | | Alternative method for defining grid farm   | | no               | None     | m           |
|                        | | offset in the x direction between rows      | | "grid"           |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``y_shear``            | | Alternative method for defining grid farm   | | no               | None     | m           |
|                        | | offset in the y direction between columns   | | "grid"           |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``min_sep_dist``       | | Minimum distance between any two turbines   | | no               | 2        | RD          |
|                        | | in a random farm                            | | "random"         |          |             |
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

To import a wind farm, set the path to a .csv file containing the per 
turbine information. In the .csv file, each column specifies a turbine
property and each row is a unique turbine. At minimum, the locations for
each turbine must be specified. Here is a small two turbine example::

         x,      y
    200.00, 0.0000
    800.00, 0.0000

Additional turbine properties can be set by adding a column with a header
equal to the yaml parameter found in the "Turbine Options" section. Here
is an example of a two turbine farm with additional properties set::

      x,       y,      HH,      yaw,      RD, thickness,   axial
    0.0,  -325.0,   110.0,   0.5236,   130.0,      13.0,    0.33
    0.0,   325.0,   110.0,  -0.5236,   130.0,      13.0,    0.33

The columns can be in any order and white space is ignored. If a property
is set in both the yaml and the imported .csv, the value in the .csv will
be used and a warning will be displayed.



Turbine Options
-----------------

This section will define all the parameters for the wind farm::

    turbines: 
        type:               <str>
        HH:                 <float>
        RD:                 <float>
        thickness:          <float>
        yaw:                <float>
        axial:              <float>
        force:              <str>
        rpm:                <float>
        read_turb_data:     <str>
        blade_segments:     <int or str>
        use_local_velocity: <bool>  
        max_chord:          <float>     
        chord_factor:       <float>     
        gauss_factor:       <float>     

+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| Option                 | Description                                   | Required (for)     | Default  | Units       |
|                        |                                               |                    |          |             |
+========================+===============================================+====================+==========+=============+
| ``type``               | | Sets the type of farm. Choices:             | yes                | None     | \-          |
|                        | | "disk" - actuator disk representation using |                    |          |             |
|                        | |          the FEniCS backend                 |                    |          |             |
|                        | | "2D_disk" - actuator disk representation    |                    |          |             |
|                        | |             optimized for 2D simulations    |                    |          |             |
|                        | | "numpy_disk" - actuator disk representation |                    |          |             |
|                        | |                that uses numpy arrays       |                    |          |             |
|                        | | "line" - actuator line representation best  |                    |          |             |
|                        | |          used with the unsteady solver      |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``HH``                 | The hub height of the turbine from ground     | all                | None     | m           |
|                        |                                               |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``RD``                 | The rotor diameter                            | all                | None     | m           |
|                        |                                               |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``yaw``                | | Determines the yaw of all turbines. Yaw is  | all                | None     | degs        |
|                        | | relative to the wind inflow direction       |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``thickness``          | The effective thickness of the rotor disk     | | "disk" or disk   | None     | m           |
|                        |                                               | | variant          |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``axial``              | The axial induction factor                    | | "disk" or disk   | None     | \-          |
|                        |                                               | | variant          |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``force``              | | the radial distribution of force            | | no               | "sine"   | \-          |
|                        | | Choices: "sine", "constant"                 | | "disk"           |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``rpm``                | | sets the revolutions per minute if using    | "line"             | 10.0     | rev/min     | 
|                        | | the alm turbine method                      |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``read_turb_data``     | | Path to .csv file with chord, lift, and     | | no               | None     | \-          |
|                        | | drag coefficients                           | | "line"           |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``blade_segments``     | | number of nodes along the rotor radius      | "line"             |"computed"| \-          |
|                        | | use "computed" to automatically set         |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``use_local_velocity`` | | use the velocity at the rotor to compute    | "line"             | True     | \-          |
|                        | | alm forces (otherwise use inflow)           |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``max_chord``          | upper limit when optimizing chord             | "line"             | 1000     | m           |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``chord_factor``       | | multiplies all the chords by a constant     | "line"             | 1.0      | \-          |
|                        | | factor                                      |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+
| ``gauss_factor``       | | factor that gets multiplied by the minimum  | "line"             | 2.0      | \-          |
|                        | | mesh spacing to set the gaussian width      |                    |          |             |
+------------------------+-----------------------------------------------+--------------------+----------+-------------+

See "Wind Farm Options" for how to specify turbine properties individually for each turbine.



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
        warp_type:         <str>
        warp_strength:     <float>
        warp_percent:      <float>
        warp_height:       <float>
        farm_num:          <int>
        farm_type:         <str>
        farm_factor:       <float>
        turbine_num:       <int>
        turbine_type:      <str>
        turbine_factor:    <float>
        refine_custom:     <list list>
        refine_power_calc: <bool>

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
| ``refine_power_calc``  | | bare minimum refinement around turbines to  |
|                        | | increase power calculation accuracy         |
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

.. note::
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
of boundary conditions: inflow, no slip, no stress. By default, inflow is 
prescribed on boundary facing into the wind, no slip on the ground and 
no stress on all other faces. These options describe the inflow boundary
velocity profile. ::

    boundary_conditions:
        vel_profile:    <str>
        HH_vel:         <float>
        vel_height:     <float, str>
        power:          <float>
        k:              <float>
        turbsim_path    <str>
        inflow_angle:   <float, list>
        boundary_names:     
            east:       <int>   
            north:      <int>   
            west:       <int>   
            south:      <int>   
            bottom:     <int>   
            top:        <int>   
            inflow:     <int>   
            outflow:    <int>   
        boundary_types:     
            inflow:     <str list> 
            no_slip:    <str list> 
            free_slip:  <str list> 
            no_stress:  <str list> 

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
| ``vel_height``         | sets the location of the reference velocity. Use "HH" for hub height                          | no           | "HH"       |
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
        type:                 <str>
        use_25d_model:        <bool>
        viscosity:            <float>
        lmax:                 <float>
        turbulence_model:     <str>
        script_iterator:      <int>             
        use_corrective_force: <bool>    
        stability_eps:        <float>             

+------------------------+--------------------------------------------------------------+--------------+---------------+
| Option                 | Description                                                  | Required     | Default       |
|                        |                                                              |              |               |
+========================+==============================================================+==============+===============+
| ``type``               | | Sets the variational form use. Choices:                    | yes          | None          |
|                        | |   "taylor_hood": Standard RANS formulation                 |              |               |
|                        | |   "stabilized": Adds a term to stabilize P1xP1 formulations|              |               |
+------------------------+--------------------------------------------------------------+--------------+---------------+
| ``viscosity``          | Kinematic Viscosity                                          | no           | 0.1           |
|                        |                                                              |              |               |
+------------------------+--------------------------------------------------------------+--------------+---------------+
| ``lmax``               | Turbulence length scale                                      | no           | 15.0          |
|                        |                                                              |              |               |
+------------------------+--------------------------------------------------------------+--------------+---------------+
| ``use_25d_model``      | | Option to enable a small amount of compressibility to mimic| | no         | False         |
|                        | | the effect of a 3D, out-of-plane flow solution in a 2D     | | "2D only"  |               |
|                        | | model.                                                     |              |               |
+------------------------+--------------------------------------------------------------+--------------+---------------+
| ``turbulence_model``   | | Sets the turbulence model.                                 | no           | mixing_length |
|                        | | Choices: mixing_length, smagorinsky, or None               |              |               |
+------------------------+--------------------------------------------------------------+--------------+---------------+
| ``script_iterator``    | debugging tool, do not use                                   | no           | 0             |
+------------------------+--------------------------------------------------------------+--------------+---------------+
|``use_corrective_force``| | add a force to the weak form to allow the inflow to recover| no           | False         |
+------------------------+--------------------------------------------------------------+--------------+---------------+
| ``stability_eps``      | | stability term to help increase the well-posedness of      | no           | 1.0           |
|                        | | the linear mixed formulation                               |              |               |
+------------------------+--------------------------------------------------------------+--------------+---------------+




Solver Options
--------------

This section lists the solver options::

    solver:
        type:              <str>
        pseudo_steady:     <bool>
        final_time:        <float>
        save_interval:     <float>
        num_wind_angles:   <int>
        endpoint:          <bool>
        velocity_path:     <str>
        power_type:        <str>
        save_power:        <bool>
        nonlinear_solver:  <str>
        newton_relaxation: <float>
        cfl_target: 0.5    <float>
        cl_iterator: 0     <int>

+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| Option                 | Description                                                    | Required (for)      | Default             |
|                        |                                                                |                     |                     |
+========================+================================================================+=====================+=====================+
| ``type``               | | Sets the solver type. Choices:                               | yes                 | None                |
|                        | |   "steady": solves for the steady state solution             |                     |                     |
|                        | |   "iterative_steady": uses iterative SIMPLE solver           |                     |                     |
|                        | |   "unsteady": solves for a time varying solution             |                     |                     |
|                        | |   "multiangle": iterates through inflow angles               |                     |                     |
|                        | |                 uses ``inflow_angle`` or [0, :math:`2\pi`]   |                     |                     |
|                        | |   "imported_inflow": runs multiple steady solves with        |                     |                     |
|                        | |                      imported list of inflow conditions      |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``pseudo_steady``      | used with unsteady solver to create a iterative steady solver. | | no                | False               |
|                        |                                                                | | "unsteady"        |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``final_time``         | The final time for an unsteady simulation                      | | no                | 1.0 s               |
|                        |                                                                | | "unsteady"        |                     |
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
| ``power_type``         | | Sets the power functional                                    | no                  | "power"             |
|                        | | Choices:                                                     |                     |                     |
|                        | |   "power": simple power calculation                          |                     |                     |
|                        | |   "2d_power": power calculation optimized for 2D runs        |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``save_power``         | | Save the power for each turbine to a text file in            | no                  | True                |
|                        | | output/``name``/data/power_data.txt                          |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``nonlinear_solver``   | | Specify the nonlinear solver type. Choices:                  | no                  | "snes"              |
|                        | |   "newton": uses the standard newton solver                  |                     |                     |
|                        | |   "snes": PETSc SNES solver                                  |                     |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``newton_relaxation``  | Set the relaxation parameter if using newton solver            | | no                | 1.0                 |
|                        |                                                                | | "newton"          |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``cfl_target``         | target CFL number for unsteady solve                           | | no                | 0.5                 |
|                        |                                                                | | "unsteady"        |                     |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+
| ``cl_iterator``        | debugging tool, do not use                                     | | no                | 0                   |
+------------------------+----------------------------------------------------------------+---------------------+---------------------+

The "multiangle" solver uses the steady solver to solve the RANS formulation.
Currently, the "multiangle" solver does not support imported domains. 


Optimization Options
--------------------

This section lists the optimization options. If you are planning on doing
optimization make sure to set ``dolfin_adjoint`` to True. ::

    optimization:
        opt_type:         <str>
        control_types:    <str list>
        layout_bounds:    <float list>
        objective_type:   <str, str list, dict>
        save_objective:   <bool>
        opt_turb_id :     <int, int list, str>
        record_time:      <str, float>
        u_avg_time:       <float>
        opt_routine:      <string>
        obj_ref:          <float>
        obj_ref0:         <float>
        taylor_test:      <bool>
        optimize:         <bool>
        gradient:         <bool>
        constraint_types: <dict>

+------------------------+----------------------------------------------------------+-----------------+--------------+
| Option                 | Description                                              | Required        | Default      |
|                        |                                                          |                 |              |
+========================+==========================================================+=================+==============+
| ``opt_type``           | Type of optimization: "minimize" or "maximize"           | no              | maximize     |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``control_types``      | | Sets the parameters to optimize. Choose Any:           | yes             | None         |
|                        | |   "yaw", "axial", "layout", "lift", "drag", "chord"    |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``layout_bounds``      | The bounding box for the layout optimization             | no              | wind_farm    |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``objective_type``     | | Sets the objective function for optimization.          | no              | power        |
|                        | | Visit :meth:`windse.objective_functions`               |                 |              |
|                        | | to see choices and additional keywords. See below to   |                 |              |
|                        | | an example for how to evaluate multiple objectives.    |                 |              |
|                        | | The first objective listed will always be used in the  |                 |              |
|                        | | optimization.                                          |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``save_objective``     | | Save the value of the objective function               | no              | True         |
|                        | | output/``name``/data/objective_data.txt                |                 |              |
|                        | | Note: power objects are saved as power_data.txt        |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``opt_turb_id``        | | Sets which turbines to optimize                        | no              | all          |
|                        | | Choices:                                               |                 |              |
|                        | |   int: optimize single turbine by ID                   |                 |              |
|                        | |   list: optimize all in list by ID                     |                 |              |
|                        | |   "all": optimize all                                  |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``record_time``        | | The amount of time to run the simulation before        | | no            | computed     |
|                        | | calculation of the objective function takes place      | | unsteady      |              |
|                        | | Choices:                                               |                 |              |
|                        | |   "computed": let the solver choose the best recording |                 |              |
|                        | |   start time based on the flow speed and domain size   |                 |              |
|                        | |   "last": only begin recording at the final_time       |                 |              |
|                        | |   <float>: time in seconds to start recording          |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``u_avg_time``         | | when to start averaging velocity for use in objective  | | no            | 5            |
|                        | | functions                                              | | unsteady      |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``opt_routine``        | | optimization method                                    | no              | SLSQP        |
|                        | | choices: SLSQP, L-BFGS-B, OM_SLSQP, SNOPT              |                 |              |
|                        | | Note: SNOPT requires custom install                    |                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``obj_ref``            | | objective reference: Sets the value of the objective   | | no            | 1.0          |
|                        | | function that will be treated as 1 by the SNOPT driver | | SLSQP         |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+
| ``obj_ref0``           | | objective reference: Sets the value of the objective   | | no            | 0.0          |
|                        | | function that will be treated as 0 by the SNOPT driver | | SLSQP         |              |
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
| ``constraint_types``   | | Allows the user to define multiple constraints.        | no              | min_dist     |
|                        | | By default, a minimum distance constraint is applied   |                 |              |
|                        | | only when performing at least a layout optimization.   |                 |              |
|                        | | additional constraints can be added similar to the way |                 |              |
|                        | | ``objective_type`` is defined. Additional detail below.|                 |              |
+------------------------+----------------------------------------------------------+-----------------+--------------+

The ``objective_type`` can be defined in three ways. First as a single string such as::

    optimization:
        objective_type: alm_power 

If the object chosen in this way has any keyword arguments, the defaults will automatically chosen. The second way is as a list of strings like::


    optimization:
        objective_type: ["alm_power", "KE_entrainment", "wake_center"]

Again, the default keyword argument will be used with this method. The final way is as a full dictionary, which allow for setting keyword arguments::

    optimization:
        objective_type:
            power: {}
            point_blockage:
                location: [0.0,0.0,240.0]
            plane_blockage_#1:
                axis: 2
                thickness: 130
                center: 240.0
            plane_blockage_#2:
                axis: 0
                thickness: 130
                center: -320.0
            cyld_kernel: 
                type: above
            mean_point_blockage:
                z_value: 240

Notice that since the objective named "power" does not have keyword arguments, an empty dictionary must be passed. For a full list of objective function visit: :meth:`windse.objective_functions`. Notice that we can have multiple version of the same objective by appending the name with "_#" and then a number. This allows us to evaluate objectives of the same type with different keyword arguments. Regardless of the number of objective types listed, currently, only the first one will be used for an optimization. 

The ``constraint_types`` option is defined in a similar way. By default the minimum distance between turbines is setup::

    constraint_types:
        min_dist:       
            target: 2   
            scale:  1   

This constraint will only be used if the ``control_types`` contains "layout". Additional constraints can be added using the same objective functions from :meth:`windse.objective_functions` by setting::

    constraint_types:
        min_dist:       
            target: 2   
            scale:  1 
        plane_blockage:
            target: 8.0
            scale: -1
            kwargs:
                axis: 2
                thickness: 130
                center: 240.0

This will still enforce the layout constraint but will additionally enforce a "plane_blockage" type constraint. By default, the constrains are setup like:

.. math::

    s * \left( c(m)-t \right) \geq 0

where :math:`c` is the constraint function, :math:`t` is the target, :math:`s` is the scale, and :math:`m` are the controls. In this configuration, we are enforcing that the result of the constraint function is greater than or equal to the target. However, we can set the scale to -1 to flip the inequality. Just like the ``objective_type``, multiple constraints of the same type can be use by appending "_#" followed by a number to the end of the name with the exception of the "min_dist" type. 