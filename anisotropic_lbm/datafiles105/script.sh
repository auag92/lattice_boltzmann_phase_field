#!/bin/bash
# print

# echo "Hello World!"
arg1=$1
# FILES=$(ls Cu_[1-5]000.dat)
# echo $FILES

gnuplot -persist <<PLOT

#set multiplot
set pm3d
unset surface
set view 0,0
set view equal xy
#set xrange [0:99]
#set yrange [0:99]
set format cb "%3.1f"
#set cbrange [0:1]
#My old color style
set palette defined (-0.00000 "red", 0.5 "green", 1 "blue")
#Gray scale color style
# set palette rgbformulae 30,31,32
#Matlab color style
#set palette rgbformulae 21,22,23
#splot 'Al_100000.dat' us 1:2:3
#replot
#set palette rgbformulae 30,31,32
splot '$arg1' us 1:2:3
# splot '$FILES' us 1:2:3
# set palette rgbformulae 21,22,23
# splot 'beta_2000.dat' us 1:2:3

# create an image file
# set terminal postscript eps enhanced color "Helvetica, 20"
# set xlabel font "Helvetica,24,bold"
# set ylabel font "Helvetica,24,bold"
# set title "Theoretical predictions of Undercooling vs Spacing"
# set xtics font "Helvetica,24,bold"
# set ytics font "Helvetica,24,bold"
# set output "test_script1.eps"
# replot

# remove all customization:
# reset
# set terminal x11

#unset multiplot
PLOT
