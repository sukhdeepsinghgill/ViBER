#!/bin/bash

make clean
make visfitsgrid

# Directory for output data
output_dir="data"
input_dir="../1_ngrf/output_data"

# Make sure the output directory exists
mkdir -p $output_dir

# Loop for 50 realizations
for realization_number in {1..100}
do
    #./visfitsgrid ../1_ngrf/output_data_pk2/Uimage_Q.beam.In_1c_${realization_number}.FITS data/Usingle_chan_${realization_number}.fits 0
    ./visfitsgrid $input_dir/image_Q.beam.In_1c_${realization_number}.FITS $output_dir/single_chan_${realization_number}.fits 0
    echo "${realization_number} Done"
	
done

echo "All realizations completed."



