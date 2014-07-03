set term pdfcairo enh color dashed font "CMR10,20"

set output "concentration.pdf"

set style line 1 lt 1 lw 3 lc rgb "#000000"
set style line 2 lt 1 lw 3 lc rgb "#006efe"
set style line 3 lt 1 lw 3 lc rgb "#c40000"

set style line 11 lt 2 lw 3 lc rgb "#000000"
set style line 12 lt 2 lw 3 lc rgb "#006efe"
set style line 13 lt 2 lw 3 lc rgb "#c40000"


set border lw 2

set ytics 0.1

set xlabel 'time {/CMMI10 t}'
set ylabel 'concentration {/CMMI10 c}'

set key  bottom right

LABEL = 'A + B {/CMSY10 \044} C'
set obj 10 rect at 7,0.57 size char strlen(LABEL), char 1 
set obj 10 front fc rgb "#FFFFFF"
set label 10 at 7,0.57 LABEL front center


set multiplot

plot [0:10][0.1:.6]\
     'without' u 1:2 w l ls 1 t 'A',\
     'without' u 1:3 w l ls 2 t 'B',\
     'without' u 1:4 w l ls 3 t 'C',\
     'with'    u 1:2 w l ls 11 not,\
     'with'    u 1:3 w l ls 12 not,\
     'with'    u 1:4 w l ls 13 not

set key bottom left height 0.7
plot [0:10][0.1:.6]\
     'without' u 1:2 w l ls 1 t 'no pressure gradient',\
     'with'    u 1:2 w l ls 11 t 'with pressure gradient'

unset multiplot