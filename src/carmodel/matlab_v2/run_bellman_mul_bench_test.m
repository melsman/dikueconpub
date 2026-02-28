% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

function [ev,t] = man_time_bellman(n,c,abar,type,pri, fipri, futhark_u0)
    if nargin < 4, type = "all"; end
    if nargin < 5, pri  = false; end
    if nargin < 6, futhark_u0 = true; end
    tstart = tic;

    ev = run_bellman_n(n,c,abar,type,pri, futhark_u0);
    
    t = toc(tstart);
    if fipri
        fprintf('Execution time: %.4f seconds\n', t);
    end
end

function t = time_bellman(n,c,abar,type,pri, fipri)

    if nargin < 4, type = "all"; end
    if nargin < 5, pri  = false; end
    if nargin < 6, fipri = false; end

    f = @() run_bellman_n(n,c,abar,type,pri);

    t = timeit(f);
    
    if fipri
        fprintf('Average execution time: %.4f seconds\n', t);
    end
end

function write_benchmark(filename, n, c, abar, method, runtime)

    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end

    % Write formatted line
    fprintf(fid, '%d %d %d %s %.6f\n', n, c, abar, method, runtime);

    % Close file
    fclose(fid);
end

function write_ev_max(filename, n, c, abar, method, ev)
    % Open file in append mode
    fid = fopen(filename, 'a');

    if fid == -1
        error('Could not open file.');
    end
    m = max(ev, [], 1);
    fprintf(fid, '# n=%d c=%d abar=%d method=%s\n', n, c, abar, char(method));
    fprintf(fid, '%.6f\n', m);
    fprintf(fid, '\n');

    fclose(fid);
end

n = 4;
abar = 20;

fid = fopen('..\fut\matlab.dat', 'w');
fclose(fid);
fid2 = fopen('..\fut\matlab_res.dat', 'w');
fclose(fid2);

for c = 5:5:35
    [ev, t] = man_time_bellman(n,c,abar,"sa",false, false, true);
    write_benchmark('..\fut\matlab.dat', n, c, abar, 'B', t);
    write_ev_max('..\fut\matlab_res.dat', n, c, abar, 'B', ev);
end

n = 4;
abar = 25;

fid = fopen('..\fut\matlabj.dat', 'w');
fclose(fid);
fid2 = fopen('..\fut\matlabj_res.dat', 'w');
fclose(fid2);

for c = 5:5:35
    [ev, t] = man_time_bellman(n,c,abar,"poly",false, false, true);
    write_benchmark('..\fut\matlabj.dat', n, c, abar, 'BJ', t);
    write_ev_max('..\fut\matlabj_res.dat', n, c, abar, 'B', ev);
end

abar = 25;
c = 35;

fid = fopen('..\fut\matlab2.dat', 'w');
fclose(fid);
fid2 = fopen('..\fut\matlab2_res.dat', 'w');
fclose(fid2);
for n = 2:2:10
    [ev, t] = man_time_bellman(n,c,abar,"poly",false, false, true);
    write_benchmark('..\fut\matlab2.dat', n, c, abar, 'B', t);
    write_ev_max('..\fut\matlab2_res.dat', n, c, abar, 'B', ev);
end


% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

