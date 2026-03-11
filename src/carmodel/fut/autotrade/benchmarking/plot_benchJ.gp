reset
set xlabel "max age"
set ylabel "time (s)"
set title "BellmanJ benchmark (n=2, c=20, N=100)"
set term pdf
set output "bench_bellmanJ.pdf"

plot "bench_bellmanJ.dat" u 4:($2==2 && $3==20 && $5==100 ? (stringcolumn(6) eq "C" ? $7 : 1/0) : 1/0) title "Futhark c" with lines, \
     ""                  u 4:($2==2 && $3==20 && $5==100 ? (stringcolumn(6) eq "M" ? $7 : 1/0) : 1/0) title "Futhark multicore" with lines, \
     ""                  u 4:($2==2 && $3==20 && $5==100 ? (stringcolumn(6) eq "O" ? $7 : 1/0) : 1/0) title "Futhark opencl" with lines, \
     "matlab_benchJ.dat" u 4:($2==2 && $3==20 && $5==100 ? (stringcolumn(6) eq "B" ? $7 : 1/0) : 1/0) title "Matlab" with lines