addpath('../matlabinclude');
addpath('../autotrade');

out_dir = '../../fut2/matlab_results_for_validation';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

n = 2;
abar = 25;
transcost = 0;
acc_0 = -5;

for c = 1:2:13
    [sol, mp, s] = bench_run_illustrations(n, c, abar, acc_0, transcost);
    write_validation(out_dir, n, c, abar, acc_0, transcost, sol, mp, s);
end

for c = 17:4:21
    [sol, mp, s] = bench_run_illustrations(n, c, abar, acc_0, transcost);
    write_validation(out_dir, n, c, abar, acc_0, transcost, sol, mp, s);
end

c = 7;

for n = 2:1:6
    [sol, mp, s] = bench_run_illustrations(n, c, abar, acc_0, transcost);
    write_validation(out_dir, n, c, abar, acc_0, transcost, sol, mp, s);
end

n = 2;

for abar = 10:5:35
    [sol, mp, s] = bench_run_illustrations(n, c, abar, acc_0, transcost);
    write_validation(out_dir, n, c, abar, acc_0, transcost, sol, mp, s);
end

function write_validation(out_dir, n, c, abar, acc_0, transcost, sol, mp, s)
    % Build [c][Ax] matrix: age 0 = mp.pnew{j}, ages 1..abar-1 = sol.p(s.ip{j}).
    prices = zeros(c, abar);
    for j = 1:c
        prices(j, 1) = mp.pnew{j};
        prices(j, 2:abar) = sol.p(s.ip{j})';
    end

    % Filename mirrors validate_eqb_all (validate_solve-N-C-A-AC-T); the
    % positive acc0 in the name matches Futhark's EQB_DATA encoding.
    acc0_pos = -acc_0;
    fname = fullfile(out_dir, sprintf('matlab-validate_solve-%d-%d-%d-%d-%d.dat', ...
                                      n, c, abar, acc0_pos, transcost));
    writematrix(prices, fname);
end
