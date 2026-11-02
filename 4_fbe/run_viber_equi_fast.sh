#!/bin/bash

make clean

make viber_equi_fast 


# Directory for output data
#output_dir="output_data"

# Make sure the output directory exists
#mkdir -p $output_dir

# Loop for 50 realizations
for realization_number in {1..500}
do
    make clean
    make viber_equi_fast 
    ./viber_equi_fast ${realization_number} GV_data_${realization_number}
    echo "GV ${realization_number} done"
done



echo "All realizations completed."


