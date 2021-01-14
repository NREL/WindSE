WindSE: Wind Systems Engineering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simple Description:
===================

WindSE is a python package that uses a FEniCS backend to perform wind farm simulations and optimization. Documentation can be found at: https://windse.readthedocs.io/en/latest/ 

Installation
============

It is easies to run WindSE within a conda environment. To install conda check this link: `Conda Installation. <https://conda.io/projects/conda/en/latest/user-guide/install/>`_ Additionally, WindSE has been tested on MacOS Catalina (10.15), but in theory should also run on linux. Windows is not recommended. 

Source Conda Installation (Script):
-----------------------------------

The easiest way to install windse is to run::

    sh install.sh <enviroment_name>

Then the enviroment can be activated using::

    conda activate <enviroment_name>

Source Conda Installation (Manual):
-----------------------------------

If you want to use the latest version or just want to setup the environment manually, follow these steps. After conda is installed, create a new environment using::

    conda create --name <enviroment_name>

You can replace the name <enviroment_name> with a different name for the environment if you want. Next we activate the environment using::

    conda activate <enviroment_name>

or whatever you named your environment. Now we need to install the dependent packages using::

    conda install -c conda-forge fenics=2019.1.0 dolfin-adjoint mshr matplotlib scipy pyyaml memory_profiler pytest

Next, we need to install the `tsfc form compilers: <https://fenics.readthedocs.io/projects/ffc/en/latest/installation.html>`_::

    pip install git+https://github.com/blechta/tsfc.git@2018.1.0
    pip install git+https://github.com/blechta/COFFEE.git@2018.1.0
    pip install git+https://github.com/blechta/FInAT.git@2018.1.0
    pip install singledispatch networkx pulp

Finally, download/clone the WindSE repo and run::

    pip install -e .

in the root folder. 

Conda-Forge Installation (Automatic):
-------------------------------

The package is available on conda-forge. To install conda check out this link: `Conda Installation. <https://conda.io/projects/conda/en/latest/user-guide/install/>`_ After conda is installed, you can automatically setup the WindSE environment using::

    conda create --name <enviroment_name>
    conda activate <enviroment_name>
    conda install -c conda-forge windse


