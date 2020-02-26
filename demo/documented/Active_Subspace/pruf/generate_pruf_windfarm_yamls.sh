#!/bin/bash

# The following file generates the YAML file with the correct wind farm information
# file. Not a fan of this but this is the quickest way to do this.

nsamples=1000
distribution_type="uniform" # "normal"
output_directory=./yaml_files/
loop_upper_bound=$(($nsamples-1))
echo $loop_upper_bound
for i in $(seq 0 $loop_upper_bound)
do
  txt_fname="windfarm_pruf_${distribution_type}_$i.txt" # construct filename for the windfarm txt
  yaml_fname="${output_directory}2d_wind_farm_PRUF_${distribution_type}_$i.yaml" # construct the yaml fname
  cp 2d_wind_farm_PRUF.yaml $yaml_fname # duplicate the pristine yaml file
  sed -i "" "11 s/pruf_wind_farm.txt/$txt_fname/" "$yaml_fname"
done
