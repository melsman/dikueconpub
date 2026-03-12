
## Benchmarking Matlab and Futhark Performance for Dynamic Programming Problem - Run_Bellman_Mul_Bench_Test

The Matlab file Run_Bellman_Mul_Branch_Test benchmarks the time and results for solving the Dynamic Programming problem in the inner loop
for a range of numbers of consumer types (n), cartypes (c) and maximum car ages (abar).

Both solving through Successive Approximation and Newton-Kantorovich-iterations are benchmarked. The benchmarked results are directly
saved in the "fut" folder after running "run_bellman_mul_bench_test". The benchmark times of the MatLab results can be compared with futhark runtimes by using "make plots.pdf" in the "fut" folder.

Currently, the C-compiled Futhark performance is equivalent to the MatLab performance for Successive Approximation. However, the Newton-Kantorovich iterations are significantly slower in Futhark.

## Benchmarking Matlab and Futhark Performance for calculating Bellman and durivative - bench_bellmanJ_test

The Matlab file bench_bellmanJ_test benchmarkes the time and results for solving 100 bellman iterations (including calculating the derivative) for 20 cartypes and a varying number of maximum car age.

Benchmarked results are saved in the "fut/autotrade/benchmarking" folder after running "bench_bellmanJ_test". Futhark benchmarks can be gotten by running "make bench_bellmanJ_all" in "fut/autotrade/benchmarking". The plot of the benchmark times can then be created by using "make bench_bellmanJ_plot" in "fut/autotrade".

Currently, the C-compiled Futhark performance is around equivalent to the MatLab performance. Therefore, I suspect the reason for the worse performancr for the Futhark in NK-iterations either has to do with some parameters that I have forgotten to set equal to each other when benchmarking, or something to do with the dpsolver code.

## Benchmarking Matlab and Futhark Performance for Newton step - bench_Newton_step_test

The Matlab file "bench_Newton_step_test" benchmarks the time for solving 60 Newton-method steps for 1 household and a maximum car age of 20 for a varying number of car types.

Benchmarked results are saved in the "fut/autotrade/benchmarking" folder after running "bench_Newton_step_test". Futhark benchmarks can be gotten by running "make bench_bellmanJ_all" in "fut/autotrade/benchmarking". The plot of the benchmark times can then be created by using "make bench_newton_single_plot" in "fut/autotrade/benchmarking".

Currently, the C-compiled Futhark performance is worse than the MatLab performance, as is the Multicore Futhark performance. Therefore, I suspect that something in the dp.solver in the Newton method is the reason for the worse Futhark performance for NK-iterations.

## Benchmarking Matlab and Futhark Performance for Newton step up to linear solve - bench_solve_linear_test

The Matlab file "bench_solve_linear_test" benchmarks the time for performing a single Newton step up to the point of solving the linear equation system.

Benchmarked results are saved in the "fut/autotrade/benchmarking" folder after running "bench_solve_linear_test". Futhark benchmarks can be gotten by running "make bench_bellmanJ_all" in "fut/autotrade/benchmarking". The plot of the benchmark times can then be created by using "make bench_linear_solver_plot" in "fut/autotrade/benchmarking".

Currently, the C-compiled Futhark performance is worse than the MatLab performance, as is the Multicore Futhark performance. Therefore, I the issue is potentially just the "lu.ols" Futhark solver being worse than the Matlab linear system solver.