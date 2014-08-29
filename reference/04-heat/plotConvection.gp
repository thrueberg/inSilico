set term pngcairo enhanced color transparent size 480, 320
set output 'uVsTime.png'

#------------------------------------------------------------------------------
# line styles
set style line 1 lt 1 pt  9 lw 3 ps 0.75 lc rgb "#000000"
set style line 2 lt 1 pt 13 lw 3 ps 0.75 lc rgb "#006efe"
set style line 3 lt 1 pt 18 lw 3 ps 0.75 lc rgb "#c40000"
set style line 4 lt 1 pt 20 lw 3 ps 0.75 lc rgb "#004989"

set style line 6 lt 1       lw 1         lc rgb "#000000"


set xlabel 'time {/CMMI10 t}'
set ylabel 'monitored values'
set key bottom right

set title "Solution at (0.5,0.5)"

plot 'convectionOutRef' using 1:2 w l  ls 2 t '{/CMMI10 u}',\
     'convectionOutRef' using 1:3 w l  ls 3 t '{/CMMI10 u}_{,{/CMMI10 x}}',\
     'convectionOutRef' using 1:4 w l  ls 4 t '{/CMMI10 u}_{,{/CMMI10 y}}'