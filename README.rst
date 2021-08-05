WindSE: Wind Systems Engineering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simple Description:
===================

WindSE is a python package that uses a FEniCS backend to perform wind farm simulations and optimization. Documentation can be found at: https://windse.readthedocs.io/en/latest/ 

Quick Start-Up Guide:
=====================

It is easiest to run WindSE within a conda environment. To install conda check this link: `Conda Installation <https://conda.io/projects/conda/en/latest/user-guide/install/>`_. Additionally, WindSE has been tested on MacOS (Catalina 10.15) and Linux (CentOS 7). Windows is not recommended. 

Quick Installation Instructions:
--------------------------------

The easiest way to install windse is to run::

    sh install.sh <enviroment_name>

Then the enviroment can be activated using::

    conda activate <enviroment_name>

Quick Demo Instructions:
------------------------

Activate the conda environment using::

    conda activate <enviroment_name>

Then to run a simple demo, navigate to <windse root>/demos/documented/Yaml_Examples/ and run::

    windse run 0-wind_farm_2D.yaml

The output of this simulation will be located in the output/2_5D_Wind_Farm/ folder. Use `Paraview <https://www.paraview.org/download/>`_ to visualize the results in the solutions/ folder. To learn what parameter can be set in the yaml file, head to the `Parameter Documentation <https://windse.readthedocs.io/en/latest/params.html>`_.



