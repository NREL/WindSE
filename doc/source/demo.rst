Demos
=====


Example Parameter Files
-----------------------

These examples show how to use the parameters file. See :doc:`params` page for more details.
All of these examples can be run using ``windse_driver run <file>``. Some file require
inputs, which can be downloaded from the section header.

1. 2D Steady State Solves (Investigates different wind farm layouts)
    a. :download:`Grid Wind Farm <demos/2D_Steady/grid_farm.yaml>`
    b. :download:`Random Wind Farm <demos/2D_Steady/random_farm.yaml>`
    c. :download:`Imported Wind Farm <demos/2D_Steady/imported_farm.yaml>`
2. 3D Steady State Solves (Investigates different domain types) (:download:`Input Data <demos/3D_Steady/Input_Data.zip>`)
    a. :download:`Box Domain <demos/3D_Steady/box.yaml>`
    b. :download:`Cylinder Domain <demos/3D_Steady/cylinder.yaml>`
    c. :download:`Imported Domain <demos/3D_Steady/imported.yaml>`
    d. :download:`Interpolated Domain <demos/3D_Steady/interpolated.yaml>`
3. Multi-Angle Solver (:download:`Input Data <demos/3D_Multi-Angle/Input_Data.zip>`)
    a. :download:`Gaussian Hill <demos/3D_Multi-Angle/params.yaml>`
4. Optimization (:download:`Input Data <demos/Optimization/Input_Data.zip>`)
    a. :download:`Axial  <demos/Optimization/axial.yaml>`
    b. :download:`Yaw <demos/Optimization/yaw.yaml>`
    c. :download:`Layout <demos/Optimization/layout.yaml>`
    d. :download:`All <demos/Optimization/all.yaml>`
5. MultiAngle Optimization (:download:`Input Data <demos/MultiAngle_Optimization/Input_Data.zip>`)
    a. :download:`2 Turbine Hill  <demos/MultiAngle_Optimization/2_turbine_hill.yaml>`
    b. :download:`3x3 Skew Hill <demos/MultiAngle_Optimization/3x3_skew.yaml>`


Example Driver Files
--------------------

These examples show how you build a custom driver if desired. Check the :doc:`api`
for details on the available functions.

1. Constructing a Gridded Wind Farm on a 2D rectangular domain: :doc:`2D_Grid
   <demos/2D_Grid_with_driver/2D_Grid.py>`.


All Demos
---------
.. toctree::
   :maxdepth: 1

   demos/2D_Grid_with_driver/2D_Grid.py.rst

