reset
set xlabel "car types"
set ylabel "time (s)"
set title "Equilibrium solve (2 households, maxage 25 years)"
set term pdf
set output "plot_eqb_cars_long.pdf"
plot "bench_solve_cars_long.dat" u 2:(stringcolumn(6) eq "C"?$7:1/0) title "Futhark c" with lines, \
     "" u 2:(stringcolumn(6) eq "M"?$7:1/0) title "Futhark multicore" with lines, \
     "" u 2:(stringcolumn(6) eq "O"?$7:1/0) title "Futhark opencl" with lines, \
     "" u 2:(stringcolumn(6) eq "U"?$7:1/0) title "Futhark cuda" with lines, \
     "matlab_eqb_cars_long.dat" u 2:6 title "Matlab" with lines
