#!/bin/bash
# get all local coordinates

folder=worker_coords
if [ -d $folder ]
then
    rm -r worker_coords
    echo "removed $folder"
fi
mkdir worker_coords
cat worker???/03_local_coords.xyz > worker_coords/local_coords.xyz

echo "All local coordinates are in folder" $folder




