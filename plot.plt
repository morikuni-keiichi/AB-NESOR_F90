reset
set format y "%.0e"
set xlabel "Number of outer iterations"
set ylabel "2-norm of relative residual"
set logscale y
set grid
set multiplot
set key spacing 0 height 1
plot "reshis.dat" title "ABGMRES-NESOR" w l
set nomultiplot
pause -1
