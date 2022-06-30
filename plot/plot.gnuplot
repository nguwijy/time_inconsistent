# Here we use argument for convenience
# this can be called from terminal using i.e.
# gnuplot -c plot.gnuplot 3
# the numbering here refers to the Figure x in
# the numerical paper

set term pdfcairo crop enhanced
set autoscale fix
set colors classic

# figure 1a
set output "gnuplot/figure1a.pdf"
set ticslevel -0.0
set key at 0,0.08,105
set xtics .01
set xtics offset graph 0,-0.05
set ytics .02
set ytics offset graph 0.05,0
set ztics 10
set xrange [-0.02:0.02]
set zrange [0:100]
set xlabel "Control" font ",11"
set ylabel "Time" font ",11"
set zlabel "Density" font ",11"
set palette rgb 33,13,10
set zlabel rotate
unset colorbox
set view 62,56
set grid
splot 'quadratic_exact_sol.txt' using 1:2:($2==0?0:$2>0.0983?0:$3>145?1/0:$3>100?100:$3) notitle with pm3d
set key at 0,0.1,95

# figure 1b
reset
set output "gnuplot/figure1b.pdf"
set ticslevel -0.0
set key at 0,0.08,105
set xtics .01
set xtics offset graph 0,-0.05
set ytics .02
set ytics offset graph 0.05,0
set ztics 10
set xrange [-0.02:0.02]
set zrange [0:100]
set xlabel "Control" font ",11"
set ylabel "Time" font ",11"
set zlabel "Density" font ",11"
set palette rgb 33,13,10
set zlabel rotate
unset colorbox
set view 62,56
splot 'quadratic_20person.txt' using 1:2:($2<0.000005?0:$2>0.0983?1/0:$3>=100?100:$3) notitle with pm3d
set key at 0,0.1,95

# figure 2
reset
set output "gnuplot/figure2.pdf"
set multiplot layout 2,2
set xlabel "Control"
set ylabel "CDF"
set key left top
set key font ",6"
set key Left
set key width -4
set style circle radius graph 0.01
set ytics 0.2
set xtics nomirror
set ytics nomirror
set grid
set xrange [-0.15:0.15]
plot 'quadratic_cdf.txt' using 1:($2==0.02?$7:1/0) with circles linecolor rgb 'black' title "Solution of 20-person; t = 0.02", 'quadratic_cdf.txt' using 1:($2==0.02?$9:1/0) with lines linecolor rgb 'blue' title "Equilibrium control; t = 0.02"
set xrange [-0.11:0.11]
set size 0.5, 0.5
plot 'quadratic_cdf.txt' using 1:($2==0.04?$7:1/0) with circles linecolor rgb 'black' title "Solution of 20-person; t = 0.04", 'quadratic_cdf.txt' using 1:($2==0.04?$9:1/0) with lines linecolor rgb 'blue' title "Equilibrium control; t = 0.04"
set xrange [-0.07:0.07]
plot 'quadratic_cdf.txt' using 1:($2==0.06?$7:1/0) with circles linecolor rgb 'black' title "Solution of 20-person; t = 0.06", 'quadratic_cdf.txt' using 1:($2==0.06?$9:1/0) with lines linecolor rgb 'blue' title "Equilibrium control; t = 0.06"
set xrange [-0.038:0.038]
set size 0.5, 0.5
plot 'quadratic_cdf.txt' using 1:($2==0.08?$7:1/0) with circles linecolor rgb 'black' title "Solution of 20-person; t = 0.08", 'quadratic_cdf.txt' using 1:($2==0.08?$9:1/0) with lines linecolor rgb 'blue' title "Equilibrium control; t = 0.08"
unset multiplot

# figure 3
reset
set output "gnuplot/figure3.pdf"
set multiplot layout 2,2
set xlabel "Control"
set ylabel "CDF"
set key font ",6"
set key Left
set key width -4
set key left top
set style circle radius graph 0.01
set xtics nomirror
set ytics nomirror
set ytics 0.2
set grid
plot 'quadratic_cdf.txt' using 1:($2==0.02?$6:1/0) with circles linecolor rgb "black" title "n = 20; m = 20; t = 0.02", 'quadratic_cdf.txt' using 1:($2==0.02?$7:1/0) with lines linecolor rgb 'blue' title "Solution of 20-person; t = 0.02"
set xrange [-0.11:0.11]
set size 0.5, 0.5
plot 'quadratic_cdf.txt' using 1:($2==0.04?$6:1/0) with circles linecolor rgb "black" title "n = 20; m = 20; t = 0.04", 'quadratic_cdf.txt' using 1:($2==0.04?$7:1/0) with lines linecolor rgb 'blue' title "Solution of 20-person; t = 0.04"
set xrange [-0.07:0.07]
plot 'quadratic_cdf.txt' using 1:($2==0.06?$6:1/0) with circles linecolor rgb "black" title "n = 20; m = 20; t = 0.06", 'quadratic_cdf.txt' using 1:($2==0.06?$7:1/0) with lines linecolor rgb 'blue' title "Solution of 20-person; t = 0.06"
set xrange [-0.038:0.038]
set size 0.5, 0.5
plot 'quadratic_cdf.txt' using 1:($2==0.08?$6:1/0) with circles linecolor rgb "black" title "n = 20; m = 20; t = 0.08", 'quadratic_cdf.txt' using 1:($2==0.08?$7:1/0) with lines linecolor rgb 'blue' title "Solution of 20-person; t = 0.08"
unset multiplot

# figure 4
reset
set output "gnuplot/figure4.pdf"
set xtics 1,2,20
set xtics nomirror
set xtics font ",15"
set yrange [0.054:0.05505]
set ytics font ",15"
set ytics nomirror
set xlabel "m" font ",19"
set ylabel "Value function" font ",19"
set key right bottom
set key font ",6"
set key Left
set key samplen 10
set key spacing 1.5 font ",14"
set grid
plot 'quadratic_value_func.txt' using 2:($1==20?$4:1/0) with lines title "n = 20" lw 2 dt (30, 6, 2, 6), 'quadratic_value_func.txt' using 2:($1==15?$4:1/0) with lines title "n = 15" lw 2 dt 3 linecolor "black", 'quadratic_value_func.txt' using 2:($1==10?$4:1/0) with lines title "n = 10" lw 2 dt (10, 6), 'quadratic_value_func.txt' using 2:($1==5?$4:1/0) with lines title "n = 5" lw 2 dt 1

# figure 5
reset
set output "gnuplot/figure5.pdf"
set xtics 1,2,20
set xtics font ",15"
set ytics font ",15"
set xtics nomirror
set ytics nomirror
set xlabel "m" font ",19"
set key Left
set key right top
set key samplen 10
set key spacing 1.5 font ",14"
set yrange [0.2:2]
set ylabel "Relative error (%)" font ",19"
set grid
plot 'quadratic_value_func.txt' using 2:($1==5?$8*100:1/0) with lines title "n = 5" lw 2 dt (30, 6, 2, 6), 'quadratic_value_func.txt' using 2:($1==10?$8*100:1/0) with lines title "n = 10" lw 2 dt 3 linecolor "black", 'quadratic_value_func.txt' using 2:($1==15?$8*100:1/0) with lines title "n = 15" lw 2 dt (10, 6), 'quadratic_value_func.txt' using 2:($1==20?$8*100:1/0) with lines title "n = 20" lw 2 dt 1

# figure 6
reset
set output "gnuplot/figure6.pdf"
set xtics 1,2,20
set xtics font ",15"
set ytics font ",15"
set xtics nomirror
set ytics nomirror
set xlabel "m" font ",19"
set ylabel "Value function" font ",19"
set key Left
set key samplen 10
set key spacing 1.5 font ",14"
set grid
plot 'quartic_value_func.txt' using 2:($1==5?$4:1/0) with lines title "n = 5" lw 2 dt (30, 6, 2, 6), 'quartic_value_func.txt' using 2:($1==10?$4:1/0) with lines title "n = 10" lw 2 dt 3 linecolor "black", 'quartic_value_func.txt' using 2:($1==15?$4:1/0) with lines title "n = 15" lw 2 dt (10, 6), 'quartic_value_func.txt' using 2:($1==20?$4:1/0) with lines title "n = 20" lw 2 dt 1

# figure 7
reset
set output "gnuplot/figure7.pdf"
set multiplot layout 2,2
set xlabel "Control"
set ylabel "CDF"
set key font ",7"
set key Left
set key width -4
set key left top
set style circle radius graph 0.01
set yrange [-.005:1.01]
set xtics nomirror
set ytics nomirror
set grid
set xtics .03
set xrange [-0.08:0.08]
plot 'quartic_cdf.txt' using 1:($2==0.02?$6:1/0) with lines linecolor rgb "blue" title "n = 20; m = 20; t = 0.02"
set xtics .02
set xrange [-0.048:0.045]
set size 0.5, 0.5
plot 'quartic_cdf.txt' using 1:($2==0.04?$6:1/0) with lines linecolor rgb "blue" title "n = 20; m = 20; t = 0.04"
set xtics .01
set yrange [-.005:1.005]
set xrange [-0.02:0.02]
plot 'quartic_cdf.txt' using 1:($2==0.06?$6:1/0) with lines linecolor rgb "blue" title "n = 20; m = 20; t = 0.06"
set xtics .002
set xrange [-0.0045:0.0048]
set size 0.5, 0.5
plot 'quartic_cdf.txt' using 1:($2==0.08?$6:1/0) with lines linecolor rgb "blue" title "n = 20; m = 20; t = 0.08"
unset multiplot
