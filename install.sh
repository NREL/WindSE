#!/bin/bash 

### Create some names
conda_base_path=$(conda info --base)
env_name=$1

### Run this to use conda in the script
source $(conda info --base)/etc/profile.d/conda.sh

### Create the Environment with mamba
conda create -y --name ${env_name} -c conda_forge python=3.9 mamba
conda activate ${env_name}

# ### Install conda-forge dependencies 
# conda install -y -c conda-forge fenics=2019.1.0=py38_9 dolfin-adjoint matplotlib scipy=1.4.1 slepc mshr hdf5 pyyaml memory_profiler pytest pytest-cov pytest-mpi coveralls pandas

# ### Install the new tsfc compiler
# pip install git+https://github.com/blechta/tsfc.git@2018.1.0
# pip install git+https://github.com/blechta/COFFEE.git@2018.1.0
# pip install git+https://github.com/blechta/FInAT.git@2018.1.0
# pip install git+https://github.com/mdolab/pyoptsparse@v1.0
# pip install singledispatch networkx pulp openmdao fatpack

# ### Install editible version of WindSE
# pip install -e .

### Install dependencies from environment
mamba env update --prefix ${conda_base_path}/envs/${env_name} --file environment.yaml --prune 
