#set term pdfcairo enh color dashed font "CMR10,20"  #6 is default!
# set output 'convergence.pdf'
set term pngcairo enhanced color transparent size 480, 320

# line style
set style line 1 lt 1 pt  9 lw 1 ps 0.75 lc rgb "#000000"


# remove borders
unset border
unset xtics
unset ytics


#------------------------------------------------------------------------------
# volume meshes
set output 'simplexMeshes.png'

set multiplot layout 1, 2;

# 2D 
set size square
plot 'squareLT.dat' w l ls 1 not

# 3D
unset ztics
set view equal xyz
splot 'cubeLT.dat' w l ls 1 not

unset multiplot

#------------------------------------------------------------------------------
# boundary meshes
set output 'boundaryMeshes.png'

set multiplot layout 1, 2;

# 2D 
set size square
plot 'squareLT_boundary.ref.dat' w l ls 1 not

# 3D
unset ztics
set view equal xyz
splot 'cubeLT_boundary.ref.dat' w l ls 1 not

unset multiplot