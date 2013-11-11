set term pdfcairo enh color dashed font "CMR10,16"  #6 is default!


#------------------------------------------------------------------------------
# line styles
set style line 1 lt 1 pt  9 lw 3 ps 0.75 lc rgb "#000000"
set style line 2 lt 1 pt 13 lw 3 ps 0.75 lc rgb "#006efe"
set style line 3 lt 1 pt 18 lw 3 ps 0.75 lc rgb "#c40000"
set style line 4 lt 1 pt 20 lw 3 ps 0.75 lc rgb "#004989"

set style line 6 lt 1       lw 1         lc rgb "#000000"


# logarithmic scale
set log xy
set format y '10^{%T}'
set border lw 2

set grid xtics

set xlabel 'mesh width {/CMMI10 h}'
set ylabel 'error'

set key bottom right

##set label 1 '1' at first 1./4.5,  first 0.2       font "CMI10,12"
##set label 2 '1' at first 1./60,   first 0.2*1e-1  font "CMI10,12"
##
##set label 3 '2' at first 1./4.5, first 0.08      font "CMI10,12"
##set label 4 '2' at first 1./60,  first 0.08*1e-2 font "CMI10,12"

set output 'interface.pdf'

plot [0.005:0.3]\
     'interface2D.dat' using (1./$1):($3) w lp ls 2 t '2D, {/CMMI10 H}^1',\
     'interface3D.dat' using (1./$1):($3) w lp ls 4 t '3D, {/CMMI10 H}^1',\
     'interface2D.dat' using (1./$1):($2) w lp ls 1 t '2D, {/CMMI10 L}_2',\
     'interface3D.dat' using (1./$1):($2) w lp ls 3 t '3D, {/CMMI10 L}_2',\
     'interface2D.dat' using (1./$1):(2./$1)       w l ls 6 not,\
     'interface2D.dat' using (1./$1):(2./$1/$1)    w l ls 6 not

