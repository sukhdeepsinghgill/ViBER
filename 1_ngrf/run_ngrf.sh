#!/bin/bash

make clean
make ngrf_2d 

# Parameters
Nx=4096
LL=1.719
ff=0.00025
spindex=0.8
freq=154.255e6
dnu=40.00e03
chan0=1.0
nchan=1

# Directory for output data
output_dir="output_data"

# Make sure the output directory exists
mkdir -p $output_dir

# Loop for 50 realizations
for realization_number in {1..50}
do
    seed=$((100100 + 413*realization_number))
    ./ngrf_2d $Nx $LL $seed $ff $spindex $freq $dnu $chan0 $nchan $output_dir/image_T_1c_${realization_number}.FITS $output_dir/image_Q.beam.In_1c_${realization_number}.FITS ${realization_number}
done

echo "All realizations completed."


