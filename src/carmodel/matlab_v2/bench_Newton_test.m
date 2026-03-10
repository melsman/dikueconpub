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



% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

