% This solves the bellman equation for a given set of prices and parameters

addpath('..\matlabinclude');
addpath('..\autotrade');

sol = bench_run_illustrations(2, 1, 25, -5, 0);
disp(mean(sol.p));

function [t,sol] = time_bench_run(n, c, abar, acc_0, transcost)

    f = @() bench_run_illustrations(n, c, abar, acc_0, transcost);
    t = timeit(f);
    sol = bench_run_illustrations(n, c, abar, acc_0, transcost);
end

t = time_bench_run(2,1,25,-5,0);
disp(t);

function write_benchmark(filename, n, c, abar, acc_0, transcost, runtime)

    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end

    % Write formatted line
    fprintf(fid, '%d %d %d %d %d %.6f\n', n, c, abar, acc_0, transcost, runtime);

    % Close file
    fclose(fid);
end


function write_avg_price(filename, n, c, abar, acc_0, transcost, sol)
    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end

    % Write formatted line
    fprintf(fid, '%d %d %d %d %d %.6f\n', n, c, abar, acc_0, transcost, mean(sol.p));

    % Close file
    fclose(fid);
end

n=2;
abar = 25;
transcost = 0;
acc_0 = -5;

fid = fopen('..\..\fut2\matlab_eqb.dat', 'w');
fclose(fid);
fid2 = fopen('..\..\fut2\matlab_eqb_res.dat', 'w');
fclose(fid2);

for c = 1:2:13
    [t, sol] = time_bench_run(n, c, abar, acc_0, transcost);
    write_benchmark('..\..\fut2\matlab_eqb.dat', n, c, abar, acc_0, transcost, t);
    write_avg_price('..\..\fut2\matlab_eqb_res.dat', n, c, abar, acc_0, transcost, sol);
end