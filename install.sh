#!/bin/bash 

### Create some names
conda_base_path=$(conda info --base)
env_name=$1

### Run this to use conda in the script
source ${conda_base_path}/etc/profile.d/conda.sh

### Create the Environment with mamba
conda create --yes --name ${env_name} --channel conda-forge python=3.9 mamba
conda activate ${env_name}

### Install dependencies from environment
mamba env update --prefix ${conda_base_path}/envs/${env_name} --file environment.yaml --prune 
