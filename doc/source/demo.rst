Demos
=====


Example Parameter Files
-----------------------

These examples show how to use the parameters file. See :doc:`params` page for more details.
All of these examples can be run using ``windse run <file>``. Some file require
inputs, which can be downloaded :download:`here <demos/Yaml_Examples/Input_Data.zip>`.

1. :download:`2D Simulations <demos/Yaml_Examples/0-wind_farm_2D.yaml>`
2. :download:`2D Layout Optimization <demos/Yaml_Examples/1-wind_farm_2D_layout_opt.yaml>`
3. :download:`3D Simulations <demos/Yaml_Examples/2-wind_farm_3D.yaml>`
4. :download:`Multi-Angle Simulations <demos/Yaml_Examples/3-multiangle_solve.yaml>`
5. :download:`Yaw Optimization <demos/Yaml_Examples/5-yaw_optimization.yaml>`
6. :download:`Multi-Angle Optimization <demos/Yaml_Examples/4-multiangle_optimization.yaml>`
7. :download:`Actuator Line Method Single-Turbine Simulation <demos/Yaml_Examples/6-alm_turbine.yaml>`

Note: These demos are extremely coarse to lower runtime for automated testing. To get better results, play with the mesh resolution and refinements. 

Example Driver Files
--------------------

These examples show how you build a custom driver if desired. Check the :doc:`api`
for details on the available functions.

1. Constructing a Gridded Wind Farm on a 2D rectangular domain: :ref:`2D Demo <demo_2d_grid>`.


Related Pages
-------------
.. toctree::
   :maxdepth: 1

   demos/Driver_Example/2D_Grid_driver.py.rst
   demos/Driver_Example/parameters.rst

