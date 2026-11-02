#!/bin/bash

# Navigate to folder1 and compile and run 
echo "********************** Processing stage 1 : ngrf **********************"
cd 1_ngrf || { echo "Failed to change directory to folder1"; exit 1; }
make clean
bash run_ngrf.sh

# Navigate to folder2 and compile and run 
echo "********************** Processing stage 2 : Visibilities **********************"
cd 2b_visfitsgrid || { echo "Failed to change directory to folder2"; exit 1; }
make clean
bash make_copies.sh
bash run_visfitsgrid.sh
cd .. || { echo "Failed to change directory back to parent directory"; exit 1; }

# Navigate to folder3 and compile and run 
echo "********************** Processing stage 3 : Convolution & gridding **********************"
cd 3_vcg || { echo "Failed to change directory to folder3"; exit 1; }
make clean
bash run_vcg.sh
cd .. || { echo "Failed to change directory back to parent directory"; exit 1; }

echo "All processes completed."

