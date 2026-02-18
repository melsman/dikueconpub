% This solves the bellman equation for a given set of prices and parameters

addpath('matlabinclude');
addpath('autotrade');

clear
clc

ev_sa, ev_poly = run_bellman_n(1, 2, 10);


% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

