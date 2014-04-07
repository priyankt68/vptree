
set title "Time taken to create vp-tree vs No. of Points"
set ylabel "Time (sec)"
set xlabel "No. of Co-ordinates"
set xrange [ 0 : 1000000000 ]
set yrange [ 0 : 500 ]
set datafile missing "-"
set xtics nomirror rotate by -65 font ",8"
set auto x
set key noenhanced
#
# First plot using linepoints
set output 'histograms.png'
set style data linespoints
set datafile separator ","
plot "../data/timings.csv" using 1:2 ti "serial", \
"../data/search_timings.csv" using 1:2 ti "parallel"
set terminal png 
set output "output.png"
replot
pause -1
