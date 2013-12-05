#!/bin/bash
clear
cp ./geo/test2geo.* ./
./test2 input.dat
mv ./*.vtk ./vtk/
rm ./test2geo.*
printf "Done \n"
