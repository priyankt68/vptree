set title "Time taken to create vp-tree vs No. of Points"
set ylabel "Time taken"
set xlabel "No. of 2D Co-ordinates"
set datafile missing "-"
#set xtics nomirror rotate by -45 font ",8"
set auto x
set key noenhanced
#
# First plot using linepoints
set style histogram
set datafile separator ","
plot "../data/timings.csv" using 1:2
pause -1
