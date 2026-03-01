% This solves the bellman equation for a given set of prices and parameters

addpath('..\matlabinclude');
addpath('..\autotrade');

clear
clc

% set parameters
mp = setparams.default(); % default parameters used for illustration

% set number of types
mp.ntypes = 2;
mp.ncartypes=2;
mp.abar_j0 = {2,2};
mp.pnew_notax = {100,100};
mp.acc_0 = {-10,-10};
mp.u_0 = {6, 6};
mp.u_a = {-0.5, -0.5};
mp.u_a_sq = {0.0, 0.0};
mp.transcost = 0;
mp.vat = 0;
mp.cartax_lo =  0;   % registration tax (below kink, K_cartax_hi)
mp.cartax_hi =  0;   % registration tax (above kink, K_cartax_hi)
mp.tax_fuel  =  0;   % proportional fuel tax 
mp.K_cartax_hi = 0;       % mp.K_cartax_hi before mp.cartax_hi tax kicks in
mp.mum = {0.1, 0.1};
mp.fixprices = true;

% update model parameter dependencies
mp = trmodel.update_mp(mp);
% model indexing
s = trmodel.index(mp);

% initialize prices where to solve bellman
% p = equilibrium.price_init(mp,s);

function [sp] = simple_prices(mp, r,s)
    sp=nan(s.np,1);
    for j=1:mp.ncartypes;
        i=1:(mp.abar_j0{j}-1);
        sp(s.ip{j}) = mp.pnew_notax{j}.*r.^i;
    end
end
p = simple_prices(mp,0.85,s);
disp(p);

ev0=zeros(s.ns, mp.ntypes); % initial guess on expected value function

% initialize arrays	% utility function
ev_sa=nan(s.ns,mp.ntypes);	% expected value function solved by successive approximations (dpsolver.sa)
ev_poly=nan(s.ns, mp.ntypes); 	% expected value function solved by poly algorithm (dpsolver.poly)

%[ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
[sol]=equilibrium.solve_once(mp, s, p);
disp(sol.ccp_scrap_tau{1});
