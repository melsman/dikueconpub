reset
set xlabel "household types"
set ylabel "time (s)"
set title "Equilibrium solve (7 car types, maxage 25 years)"
set term pdf
set output "plot_eqb_households.pdf"
plot "bench_solve_households.dat" u 1:(stringcolumn(6) eq "C"?$7:1/0) title "Futhark c" with lines, \
     "" u 1:(stringcolumn(6) eq "M"?$7:1/0) title "Futhark multicore" with lines, \
     "" u 1:(stringcolumn(6) eq "O"?$7:1/0) title "Futhark opencl" with lines, \
     "" u 1:(stringcolumn(6) eq "U"?$7:1/0) title "Futhark cuda" with lines, \
     "matlab_eqb_households.dat" u 1:6 title "Matlab" with lines
