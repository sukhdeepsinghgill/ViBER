#!/bin/bash

make clean
make vcg

# Directory for output data
output_dir="output_data"
input_dir="../2b_visfitsgrid/data"

# Make sure the output directory exists
mkdir -p $output_dir

./vcg $input_dir/single_chan_1.fits inputparameters $output_dir/GV_data_1 1 $output_dir/sampling

for realization_number in {2..500}
do
    ./vcg $input_dir/single_chan_${realization_number}.fits inputparameters $output_dir/GV_data_${realization_number} 0 ph
    echo "############${realization_number} completed.###########################"
done

echo "All realizations completed."

##./vcg ../2b_visfitsgrid/data/Usingle_chan_${realization_number}.fits inputparameters $output_dir/UGV_data_${realization_number}

