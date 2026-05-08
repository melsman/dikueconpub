% This solves the bellman equation for a given set of prices and parameters

addpath('../matlabinclude');
addpath('../autotrade');

sol = bench_run_illustrations(2, 1, 25, -5, 0);
disp(mean(sol.p));

t = time_bench_run(2,1,25,-5,0);
disp(t);

n=2;
abar = 25;
transcost = 0;
acc_0 = -5;

cars_fid_bench = '../../fut2/matlab_eqb_cars.dat';
cars_fid_res = '../../fut2/matlab_eqb_res_cars.dat';

fid = fopen(cars_fid_bench, 'w');
fclose(fid);
fid2 = fopen(cars_fid_res, 'w');
fclose(fid2);

for c = 1:2:13
    f = @() bench_run_illustrations(n, c, abar, acc_0, transcost);
    sol = bench_run_illustrations(n, c, abar, acc_0, transcost);
    stats = man_time_it(f, 3, 0);
    write_benchmark(cars_fid_bench, n, c, abar, acc_0, transcost, stats.mean, stats.std, stats.se);
    write_avg_price(cars_fid_res, n, c, abar, acc_0, transcost, sol);
end

c = 7;

households_fid_bench = '../../fut2/matlab_eqb_households.dat';
households_fid_res = '../../fut2/matlab_eqb_res_households.dat';

fid = fopen(households_fid_bench, 'w');
fclose(fid);
fid2 = fopen(households_fid_res, 'w');
fclose(fid2);

for n = 2:1:6
    f = @() bench_run_illustrations(n, c, abar, acc_0, transcost);
    sol = bench_run_illustrations(n, c, abar, acc_0, transcost);
    stats = man_time_it(f, 3, 0);
    write_benchmark(households_fid_bench, n, c, abar, acc_0, transcost, stats.mean, stats.std, stats.se);
    write_avg_price(households_fid_res, n, c, abar, acc_0, transcost, sol);
end

age_fid_bench = '../../fut2/matlab_eqb_age.dat';
age_fid_res = '../../fut2/matlab_eqb_res_age.dat';

fid = fopen(age_fid_bench, 'w');
fclose(fid);
fid2 = fopen(age_fid_res, 'w');
fclose(fid2);

n = 2;
for abar = 10:5:35
    f = @() bench_run_illustrations(n, c, abar, acc_0, transcost);
    sol = bench_run_illustrations(n, c, abar, acc_0, transcost);
    stats = man_time_it(f, 3, 0);
    write_benchmark(age_fid_bench, n, c, abar, acc_0, transcost, stats.mean, stats.std, stats.se);
    write_avg_price(age_fid_res, n, c, abar, acc_0, transcost, sol);
end

function [t,sol] = time_bench_run(n, c, abar, acc_0, transcost)

    f = @() bench_run_illustrations(n, c, abar, acc_0, transcost);
    t = timeit(f);
    sol = bench_run_illustrations(n, c, abar, acc_0, transcost);
end

function stats = man_time_it(f, num_repeats, num_warmups)

    if nargin < 2
        num_repeats = 30;
    end

    if nargin < 3
        num_warmups = 3;
    end

    % Warm-up runs
    for k = 1:num_warmups
        f();
    end

    % Timed runs
    times = zeros(num_repeats, 1);

    for k = 1:num_repeats
        tic;
        f();
        times(k) = toc;
    end

    stats.mean = mean(times);
    stats.median = median(times);
    stats.std = std(times);                  % sample std, normalized by N-1
    stats.se = stats.std / sqrt(num_repeats);
    stats.min = min(times);
    stats.max = max(times);
    stats.num_repeats = num_repeats;
    stats.num_warmups = num_warmups;
    stats.times = times;

end

function write_benchmark(filename, n, c, abar, acc_0, transcost, runtime, std_dev, std_error)

    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end

    % Write formatted line
    fprintf(fid, '%d %d %d %d %d %.6f %.6f %.6f\n', n, c, abar, acc_0, transcost, runtime, std_dev, std_error);

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
