#!/bin/bash
./clean.sh
sleep 1
reset
cp ./geo/test2geo.* ./
./test2 input.dat
mv ./*.vtk ./vtk/
rm ./test2geo.*
printf "\nStarting Paraview..."
./view.sh
printf "Done \n"
