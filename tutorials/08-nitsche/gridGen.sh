#!/bin/bash

# -------------------------------------------------------------------------------
# Check number of input arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 (1|2|3) (S|U)"
    exit 1
fi

# Read input arguments
DIM=$1
MODE=$2

# Variables set in function of the second input argument
EXE=
PRE=
SUF=
case $MODE in
    S)
        EXE=$INSILICOROOT/tools/bin/unitCubeSGF
        PRE="grid"
        SUF=".sgf"
        ;;
    U)
        EXE=$INSILICOROOT/tools/bin/unitCubeSMFsimplex
        PRE="mesh"
        SUF=".smf"
        ;;
    *)
        echo "Mode has to have the value \"S\" or \"U\", but $2 was passed"
        exit 1
esac

# Output filename prefix
PRE+=$DIM
PRE+="D."

# ------------------------------------------------------------------------------
# Choose dimension
case $DIM in
    1) # 1D - grids
        for N in 004 006 008 012 016 024 032 048 064 096 128 192 256; do
            NAME=$PRE$N$SUF
            $EXE $N > $NAME;
         done
        ;;
    2) # 2D - grids
        for N in 004 006 008 012 016 024 032 048 064 096 128; do
            NAME=$PRE$N$SUF
            $EXE $N $N > $NAME;
        done
        ;;
    3) # 3D - grids
        for N in 004 006 008 012 016 024 032 048; do
            NAME=$PRE$N$SUF
            $EXE $N $N $N > $NAME;
        done
        ;;
    *)
        echo "The value of the input argument was $1 and has to be 1, 2, or 3"
        exit 1
 
esac
