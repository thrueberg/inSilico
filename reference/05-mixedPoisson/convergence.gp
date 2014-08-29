#set term pdfcairo enh color dashed font "CMR10,20"  #6 is default!
# set output 'convergence.pdf'
set term pngcairo enhanced color transparent size 480, 320
set output 'convergence.png'

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

set key at first 0.2,first 1.e-6  spacing 1.2

set label 1 '1' at first 1./4.5, first 0.2       font "CMI10,12"
set label 2 '1' at first 1./600,  first 0.2*1e-2  font "CMI10,12"

set label 3 '2' at first 1./4.5, first 0.08      font "CMI10,12"
set label 4 '2' at first 1./600,  first 0.08*1e-4 font "CMI10,12"

set label 5 '3' at first 1./4.5, first 0.0008      font "CMI10,12"
set label 6 '3' at first 1./600,  first 0.0008*1e-6 font "CMI10,12"

plot [0.001:0.3]\
     'convergence.dat' index 0 using (1./$1):($3) w lp  ls 2 t '{/CMMI10 p}=1, {/CMMI10 H}^1',\
     'convergence.dat' index 1 using (1./$1):3 w lp  ls 4 t '{/CMMI10 p}=2, {/CMMI10 H}_1',\
     'convergence.dat' index 0 using (1./$1):2 w lp  ls 1 t '{/CMMI10 p}=1, {/CMMI10 L}_2',\
     'convergence.dat' index 1 using (1./$1):2 w lp  ls 3 t '{/CMMI10 p}=2, {/CMMI10 L}_2',\
     'convergence.dat' index 0 using (1./$1):(1./$1)       w l ls 6 not,\
     'convergence.dat' index 0 using (1./$1):(2./$1/$1)    w l ls 6 not,\
     'convergence.dat' index 0 using (1./$1):(.1/$1/$1/$1) w l ls 6 not