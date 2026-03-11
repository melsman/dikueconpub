% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

tstart = tic;
[ev, iter] = bench_Newton(1, 20, 80, 20, 0);
t = toc(tstart);
disp(iter);
disp(max(ev));
disp(t);


function t = time_bench_Newton(n, c, abar, iter, tol)

    f = @() bench_Newton(n, c, abar, iter, tol);
    t = timeit(f);
end

function write_benchmark(filename, n, c, abar, iter, tol, runtime);

    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end

    % Write formatted line
    fprintf(fid,  '%s %d %d %d %d %d %s %.6f\n', "J", n, c, abar, iter, tol, "B", runtime);

    % Close file
    fclose(fid);
end

n = 1;
abar = 20;
tol = 0;
iter = 60;

fid = fopen('..\fut\autotrade\benchmarking\matlat_newton_bench.dat', 'w');
fclose(fid);
for c = 5:5:50
    t = time_bench_Newton(n, c, abar, iter, tol);
    write_benchmark('..\fut\autotrade\benchmarking\matlat_newton_bench.dat', n, c, abar, iter, tol, t);
end

% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

