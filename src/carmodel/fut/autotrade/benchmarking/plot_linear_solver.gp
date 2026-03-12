reset
set xlabel "car types"
set ylabel "time (s)"
set title "Linear solve benchmark (1 household, maxage 25 years)"
set term pdf
set output "bench_linear_solver.pdf"

plot "bench_linear_solver.dat" u 2:($1==1 && $3==25 ? (stringcolumn(4) eq "C" ? $5 : 1/0) : 1/0) title "Futhark c" with lines, \
     ""                        u 2:($1==1 && $3==25 ? (stringcolumn(4) eq "M" ? $5 : 1/0) : 1/0) title "Futhark multicore" with lines, \
     ""                        u 2:($1==1 && $3==25 ? (stringcolumn(4) eq "O" ? $5 : 1/0) : 1/0) title "Futhark opencl" with lines, \
     "matlab_linear_bench.dat" u 3:($2==1 && $4==25 ? (stringcolumn(5) eq "B" ? $6 : 1/0) : 1/0) title "Matlab" with lines