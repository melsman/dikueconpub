
## Benchmarking Matlab and Futhark Performance for Dynamic Programming Problem - Run_Bellman_Mul_Bench_Test

The Matlab file Run_Bellman_Mul_Branch_Test benchmarks the time and results for solving the Dynamic Programming problem in the inner loop
for a range of numbers of consumer types (n), cartypes (c) and maximum car ages (abar).

Both solving through Successive Approximation and Newton-Kantorovich-iterations are benchmarked. The benchmarked results are directly
saved in the "fut" folder after running "run_bellman_mul_bench_test". The benchmark times of the MatLab results can be compared with futhark runtimes by using "make plots.pdf" in the "fut" folder.

Currently, the C-compiled Futhark performance is equivalent to the MatLab performance for Successive Approximation. However, the Newton-Kantorovich iterations are significantly slower in Futhark.

## Benchmarking Matlab and Futhark Performance for calculating Bellman and durivative - bench_bellmanJ_test

The Matlab file bench_bellmanJ_test benchmarkes the time and results for solving 100 bellman iterations (including calculating the derivative) for 20 cartypes and a varying number of maximum car age.

Benchmarked results are saved in the "fut/autotrade" folder after running "bench_bellmanJ_test". Futhark benchmarks can be gotten by running "make bench_bellmanJ_all" in "fut/autotrade". The plot of the benchmark times can then be created by using "make bench_bellmanJ_plot" in "fut/autotrade".

Currently, the C-compiled Futhark performance is around equivalent to the MatLab performance. Therefore, I suspect the reason for the worse performancr for the Futhark in NK-iterations either has to do with some parameters that I have forgotten to set equal to each other when benchmarking, or something to do with the dpsolver code.