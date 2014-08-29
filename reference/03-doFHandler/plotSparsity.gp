set term pngcairo enhanced color transparent size 840, 320

# line style
set style line 1 lt 1 pt  9 lw 1 ps 0.75 lc rgb "#000000"


# remove borders
unset xtics
unset ytics

set yrange [] reverse

set output 'sparsity.png'

set size square

set multiplot layout 1, 4;

# mesh 
set title "mesh"
plot 'square_20.dat' w l ls 1 not

# linear
set title "linear"
plot 'sparsity.1.ref.dat' w d ls 1 not

# quadratic
set title "quadratic"
plot 'sparsity.2.ref.dat' w d ls 1 not

# cubic
set title "cubic"
plot 'sparsity.3.ref.dat' w d ls 1 not

unset multiplot

