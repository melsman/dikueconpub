% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

[ev, t] = bench_solve_linear(30, 25);
disp(t);
disp(sum(ev));
disp(max(ev));

function write_benchmark_linear(filename, n, c, abar, method, runtime)

    fid = fopen(filename, 'a');
    if fid == -1
        error('Could not open file.');
    end

    fprintf(fid, '%s %d %d %d %s %.6f\n', "J", n, c, abar, method, runtime);

    fclose(fid);
end

fid = fopen('..\fut\autotrade\benchmarking\matlab_linear_bench.dat', 'w');
fclose(fid);

n = 1;
abar = 25;

for c = 5:5:50
    [x, t] = bench_solve_linear(c, abar);
    write_benchmark_linear('..\fut\autotrade\benchmarking\matlab_linear_bench.dat', ...
                           n, c, abar, "B", t);
end

% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

