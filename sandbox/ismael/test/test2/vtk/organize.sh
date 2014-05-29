#!/bin/bash

reset

DIRECTORY="test-"

for (( N=1; c <=99; N++))
do
	#Numeration correction
	if [ "$N" -lt "10" ];then
		DIRC="$DIRECTORY""0"
	else
		DIRC="$DIRECTORY"
	fi

	#Check and Do
	if [ ! -d "$DIRC""$N" ]; then
		mkdir "$DIRC""$N"
		mv ./*.vtk "./$DIRC""$N" 2> /dev/null
		cp ../input.dat "./$DIRC""$N"

		echo Work done.
		exit
	fi
done
