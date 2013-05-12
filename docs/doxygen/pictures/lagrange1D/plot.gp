# terminal for PNG
set term pngcairo enhanced color transparent

# line styles
set style line 1 lt 1 lw 2 lc rgb "#000000"

# ------------------------------------------------------------------------------
# plot boundaries
xmin = 0.0
xmax = 1.0
ymin = -0.5
ymax = 1.1
set xrange [xmin:xmax]
set yrange [ymin:ymax]

# ------------------------------------------------------------------------------
# Lagrange functions
pMax = 3

# p = 0
l0_0(x) = ( (x < 0.) || (x > 1.) ? 1./0. : 1. )
l0(x,s) = (\
        s == 0? l0_0(x) :\
        1./0.)

# p = 1
l1_0(x) = l0_0(x) * (1. - x)
l1_1(x) = l0_0(x) * x
l1(x,s) = (\
        s == 0? l1_0(x) :\
        s == 1? l1_1(x) :\
        1./0. )

# p = 2
l2_0(x) = l0_0(x) * (1. - x) * (1. - 2.*x)
l2_1(x) = l0_0(x) * 4. * x * (1. - x)
l2_2(x) = l0_0(x) * x * (2.*x - 1.)
l2(x,s) = (\
        s == 0? l2_0(x) :\
        s == 1? l2_1(x) :\
        s == 2? l2_2(x) :\
        1./0. )

# p = 3
z0(x) = 1 - x
z1(x) = x
l3_0(x) = l0_0(x) * 0.5 * z0(x) * (3. * z0(x) - 1.) * (3. * z0(x) - 2.)
l3_1(x) = l0_0(x) * 4.5 * z0(x) * (3. * z0(x) - 1.) * z1(x);
l3_2(x) = l0_0(x) * 4.5 * z1(x) * (3. * z1(x) - 1.) * z0(x);
l3_3(x) = l0_0(x) * 0.5 * z1(x) * (3. * z1(x) - 1.) * (3. * z1(x) - 2.)
l3(x,s) = (\
        s == 0? l3_0(x) :\
        s == 1? l3_1(x) :\
        s == 2? l3_2(x) :\
        s == 3? l3_3(x) :\
        1./0. )

# ------------------------------------------------------------------------------
# generic lagrange fun in function of the degree
lag(x,p,s) = (\
           p == 0 ? l0(x,s) :\
           p == 1 ? l1(x,s) :\
           p == 2 ? l2(x,s) :\
           p == 3 ? l3(x,s) :\
           1./0. )


# title name for the key
titname( p, x ) = sprintf( '{/CMMI10 L}^%d_%d', p, s )

# key at the bottom line of the plot
set key inside bottom center horizontal samplen 1.0 width 3.0

# loop over all available degrees
do for [p=0:pMax] {
# output file in function of the degree
outfile = sprintf('lagrange_%d.png',p)
set output outfile
# a fat line along the x-axis
set arrow 1 from first 0,0 to 1,0 nohead back ls 1
# plot for all functions for this degree
plot  for [s=0:p] lag(x,p,s) t titname(p,s)
}