set title "Time taken for K-NN search (K=8)"
set ylabel "Time (sec)"
set xlabel "No. of Co-ordinates in Query Set (Q)"
set xrange [ 0 : 1000000000 ]
set yrange [ 0 : 20 ]
set datafile missing "-"
set xtics nomirror rotate by -65 font ",8"
set auto x
#set linestyle 1 lt 2 lw 3
#set key box linestyle 1
set key box left top
#
# First plot using linepoints

set style data linespoints
set datafile separator ","
plot "../data/serial_vs_parallel.csv" using 1:2 ti "parallel", "../data/serial_vs_parallel.csv" using 1:3 ti "sequential"
set terminal png 
set output "serial_vs_parallel.png"
replot
pause -1
