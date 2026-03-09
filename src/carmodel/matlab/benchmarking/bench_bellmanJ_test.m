% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

[ev, ccp, dev] = bench_bellmanJ(2,15,15,100);
%disp(ev);
%disp(dev(:,1));

function t = time_bellman(n, c, abar, N)

    f = @() bench_bellmanJ(n, c,abar,N);
    t = timeit(f);
end

t = time_bellman(2,15,30,100);
disp(t)

function write_benchmark(filename, n, c, abar, N, runtime)

    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end

    % Write formatted line
    fprintf(fid,  '%s %d %d %d %d %s %.6f\n', "J", n, c, abar, N, "B", runtime);

    % Close file
    fclose(fid);
end

fid = fopen('..\..\fut2\autotrade\matlab_benchj.dat', 'w');
fclose(fid);

n = 2;
c = 20;
N = 100;

for abar = 5:5:40
    t = time_bellman(n, c, abar, N);
    write_benchmark('..\..\fut2\autotrade\matlab_benchj.dat', n, c, abar, N, t);
end
% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

