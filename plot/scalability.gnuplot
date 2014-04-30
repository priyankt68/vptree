set title "Time taken for parallel K-NN search (K=8)"
set ylabel "Time (sec)"
set xlabel "No. of Processors (P)"
set xrange [ 1 : 8 ]
set yrange [ -1 : 5 ]
set datafile missing "-"
#set xtics nomirror rotate by -65 font ",8"
set auto x
#set linestyle 1 lt 2 lw 3
#set key box linestyle 1
#set key box left top
#
# First plot using linepoints

set style data linespoints
set datafile separator ","
plot "../data/scalability.csv" using 1:2 ti "Time taken/Processors"
set terminal png 
set output "scalability.png"
replot
pause -1
