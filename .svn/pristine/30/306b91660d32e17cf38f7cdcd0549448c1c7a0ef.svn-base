#!/bin/bash

# -------------------------------------------------------------------------------
# Check number of input arguments
if [ $# -ne 1 ]; then
    echo "Usage:  $0 (2|3) "
    exit 1
fi

# validate input
if [[ ($1 != 2) && ($1 != 3) ]]; then
    echo "Wrong 1st input parameter: $2"
    exit 1
fi

# Read input arguments
DIM=$1

# ------------------------------------------------------------------------------
# Variables set in function of the second input argument

# executable
EXE="./linearElastic"
EXE+=$DIM
EXE+="D"

# output file for convergence data
CONV="linearElastic"
CONV+=$DIM
CONV+="D"
if [[ "$EXENAME" = "box" ]] ; then CONV+=$MODE; fi
CONV+=".dat"

PRE=
SUF=".smf"
case $DIM in
    2)
        PRE="quad."
        ;;
    3)
        PRE="cube."
        ;;
    *)
        echo "Mode has to have the value \"S\" or \"U\", but $2 was passed"
        exit 1
esac


# execute in a loop
echo "#  N   L2   " > $CONV
for name in $PRE*$SUF; do
    echo $name;
    echo -e ${name:5:3} "\c" >> $CONV && $EXE $name $ARGS >> $CONV;
done
