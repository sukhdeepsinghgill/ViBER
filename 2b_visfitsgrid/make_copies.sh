#!/bin/bash

# Source file
source_file="single_chan.fits"

# Check if the source file exists
if [ ! -f "$source_file" ]; then
    echo "Source file '$source_file' does not exist!"
    exit 1
fi

# Loop to create 50 copies
for relzn in {201..500}
do
   # dest_file1="data/Usingle_chan_${relzn}.fits"
    dest_file2="data/single_chan_${relzn}.fits"
   # cp "$source_file" "$dest_file1"
    cp "$source_file" "$dest_file2"
    echo "${relzn} Copied"
done

echo "All copies created successfully."

