# Steps to perform an analysis
#
#
# 1) Adapt the file mesh/microfluidic_chip.geo (e.g. sizes, element length...)
# 2) Enter mesh directory and run gmsh
cd mesh
~trueberg/bin/gmsh -2 microfluidic_chip.geo
# 3) Convert to smf
~trueberg/src/inSilico/tools/bin/gmsh2smf microfluidic_chip.msh
# 4) Clean the file (remove orphan nodes)
~trueberg/src/inSilico/tools/bin/smfClean microfluidic_chip.triangle.smf microfluidic_chip.smf
# 5) Copy file to working folder and go back
cp microfluidic_chip.smf ../.
cd ..

