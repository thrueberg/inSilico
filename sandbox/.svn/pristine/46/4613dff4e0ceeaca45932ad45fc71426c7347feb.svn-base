#!/bin/bash

# -------------------------------------------------------------------------------
# Check number of input arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 N  D"
    exit 1
fi

# Read input arguments
NY=$1
DELTA=$2

# path to find tools
TOOLSPATH=$INSILICOROOT/tools/bin

#generate mesh
UC=unitCubeSMF
NX=$((5 * $NY / 2))
SMFTMP=unit.smf
$TOOLSPATH/$UC  $NX $NY > $SMFTMP

#perform affine map
AF=smfAffine
AFFTMP=affTmp
echo -e "5\n 0\n 0\n 0\n 2\n 0\n 0\n 0\n 1\n -1\n -1\n 0\n" > $AFFTMP

$TOOLSPATH/$AF $SMFTMP < $AFFTMP
AFFOUT=${SMFTMP%.smf}.aff.smf


#perform map
CM=./constrict
$CM $AFFOUT $DELTA

CMOUT=${AFFOUT%.smf}.map.smf

#clean
mv $CMOUT  pipe.$NY.smf
rm $SMFTMP
rm $AFFTMP
rm $AFFOUT



