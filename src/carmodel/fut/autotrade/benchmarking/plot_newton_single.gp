reset
set xlabel "car types"
set ylabel "time (s)"
set title "Time (1 household, maxage 20 years) for 60 Newton steps"
set term pdf
set output "bench_newton_single.pdf"

plot "bench_newton_single.dat" u 2:($1==1 && $3==20 && $4==60 && $5==0 ? (stringcolumn(6) eq "C" ? $7 : 1/0) : 1/0) title "Futhark c" with lines, \
     ""                        u 2:($1==1 && $3==20 && $4==60 && $5==0 ? (stringcolumn(6) eq "M" ? $7 : 1/0) : 1/0) title "Futhark multicore" with lines, \
     ""                        u 2:($1==1 && $3==20 && $4==60 && $5==0 ? (stringcolumn(6) eq "O" ? $7 : 1/0) : 1/0) title "Futhark opencl" with lines, \
     "matlat_newton_bench.dat" u 3:($2==1 && $4==20 && $5==60 && $6==0 ? (stringcolumn(7) eq "B" ? $8 : 1/0) : 1/0) title "Matlab" with lines