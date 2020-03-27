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
  eagle_directory=/projects/uqpwrsys/kpanda/windfarm_files/ # Directory where the windfarm txts are stored
  txt_fname="${eagle_directory}windfarm_pruf_${distribution_type}_$i.txt" # construct filename for the windfarm txt
  yaml_fname="${output_directory}2d_wind_farm_PRUF_${distribution_type}_$i.yaml" # construct the yaml fname
  case_name="2_5D_Wind_Farm_PRUF_${i}"
  # echo $txt_fname $case_name
  cp 2d_wind_farm_PRUF.yaml $yaml_fname # duplicate the pristine yaml file
  sed -i "" "3 s/2_5D_Wind_Farm_PRUF/${case_name}/" "$yaml_fname"
  sed -i "" "13 s~pruf_wind_farm.txt~${txt_fname}~" "$yaml_fname"

  # Uncomment the following for low-fidelity yamls
  sed -i "" "25 s/200/100/" "$yaml_fname"
  sed -i "" "26 s/200/100/" "$yaml_fname"
done

cd ./yaml_files/
tar -czf windfarm_PRUF_uniform_lofi_yamls.tar.gz ./*.yaml
# tar -czf windfarm_PRUF_uniform_yamls.tar.gz ./*.yaml
cd ../
