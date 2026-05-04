reset
set xlabel "max age"
set ylabel "time (s)"
set title "Equilibrium solve (7 car types)"
set term pdf
set output "plot_eqb_age.pdf"
plot "bench_solve_age.dat" u 3:(stringcolumn(6) eq "C"?$7:1/0) title "Futhark c" with lines, \
     "" u 3:(stringcolumn(6) eq "M"?$7:1/0) title "Futhark multicore" with lines, \
     "" u 3:(stringcolumn(6) eq "O"?$7:1/0) title "Futhark opencl" with lines, \
     "" u 3:(stringcolumn(6) eq "U"?$7:1/0) title "Futhark cuda" with lines, \
     "matlab_eqb_age.dat" u 3:6 title "Matlab" with lines
