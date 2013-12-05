# set pdf terminal with computed size
set term pdfcairo enh color dashed font "CMR10,18" 

#------------------------------------------------------------------------------
# line styles
set style line 1 lt 1 pt  9 lw 3 ps 0.75 lc rgb "#000000"
set style line 2 lt 1 pt 13 lw 3 ps 0.75 lc rgb "#006efe"
set style line 3 lt 1 pt 18 lw 3 ps 0.75 lc rgb "#c40000"
set style line 4 lt 4 pt 20 lw 3 ps 0.75 lc rgb "#004989"
set style line 5 lt 5 pt 11 lw 3 ps 0.75 lc rgb "#898900"

set output 'error.pdf'

err(x,y) = abs(x-y)

set log y

plot 'bdf1.dt0p1000' u 1:(err($2,$3)) w l ls 1 t 'BDF 1',\
     'bdf2.dt0p1000' u 1:(err($2,$3)) w l ls 2 t 'BDF 2',\
     'bdf3.dt0p1000' u 1:(err($2,$3)) w l ls 3 t 'BDF 3',\
     'am2.dt0p1000'  u 1:(err($2,$3)) w l ls 4 t 'AM  2',\
     'am3.dt0p1000'  u 1:(err($2,$3)) w l ls 5 t 'AM  3'