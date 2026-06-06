reset
set xlabel "car types"
set ylabel "time (s)"
set title "Equilibrium solve (2 households, maxage 25 years)"
set term pdf
set output "plot_eqb_cars_longest.pdf"
has_matlab = system("[ -f matlab_eqb_cars_longest.dat ] && echo 1 || echo 0") + 0
if (has_matlab) {
    plot "bench_solve_cars_longest.dat" u 2:(stringcolumn(6) eq "C"?$7:1/0) title "Futhark c" with lines, \
         "" u 2:(stringcolumn(6) eq "M"?$7:1/0) title "Futhark multicore" with lines, \
         "" u 2:(stringcolumn(6) eq "O"?$7:1/0) title "Futhark opencl" with lines, \
         "" u 2:(stringcolumn(6) eq "U"?$7:1/0) title "Futhark cuda" with lines, \
         "matlab_eqb_cars_longest.dat" u 2:6 title "Matlab" with lines
} else {
    plot "bench_solve_cars_longest.dat" u 2:(stringcolumn(6) eq "C"?$7:1/0) title "Futhark c" with lines, \
         "" u 2:(stringcolumn(6) eq "M"?$7:1/0) title "Futhark multicore" with lines, \
         "" u 2:(stringcolumn(6) eq "O"?$7:1/0) title "Futhark opencl" with lines, \
         "" u 2:(stringcolumn(6) eq "U"?$7:1/0) title "Futhark cuda" with lines
}
