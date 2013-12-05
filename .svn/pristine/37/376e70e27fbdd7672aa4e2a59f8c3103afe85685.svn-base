#!/bin/bash

# -------------------------------------------------------------------------------
# Check number of input arguments
if [ $# -ne 3 ]; then
    echo "Usage:  $0 (box|interface) (2|3) (S|U)"
    exit 1
fi

# validate input
if [[ ("$1" != "box") && ("$1" != "interface") ]]; then
    echo "Wrong 1st input parameter: $1"
    exit 1
fi

if [[ ($2 != 2) && ($2 != 3) ]]; then
    echo "Wrong 2nd input parameter: $2"
    exit 1
fi

if [[ ("$3" != "S") && ("$3" != "U") ]]; then
    echo "Wrong 3rd input parameter: $3"
    exit 1
fi

# Read input arguments
EXENAME=$1
DIM=$2
MODE=$3

# ------------------------------------------------------------------------------
# Variables set in function of the second input argument

# executable
EXE="./"
EXE+=$EXENAME$DIM
EXE+="D"
if [[ "$EXENAME" = "box" ]] ; then EXE+=$MODE; fi

# output file for convergence data
CONV=$1
CONV+=$DIM
CONV+="D"
if [[ "$EXENAME" = "box" ]] ; then CONV+=$MODE; fi
CONV+=".dat"

PRE=
SUF=
case $MODE in
    S)
        PRE="grid"
        SUF=".sgf"
        ;;
    U)
        PRE="mesh"
        SUF=".smf"
        ;;
    *)
        echo "Mode has to have the value \"S\" or \"U\", but $2 was passed"
        exit 1
esac

# command line arguments (only for interface program)
ARGS=
if [[ "$EXENAME" = "interface" ]]; then ARGS+="0.2222 3. 1."; fi

# Output filename prefix
PRE+=$DIM
PRE+="D."

echo "#  N   L2   H1   num-CG-iter " > $CONV
for name in $PRE*$SUF; do
    echo $name;
    echo -e ${name:7:3} "\c" >> $CONV && $EXE $name $ARGS >> $CONV;
done
