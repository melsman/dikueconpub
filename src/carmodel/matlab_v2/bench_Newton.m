% This solves the bellman equation for a given set of prices and parameters


function [ev, iter_n] = bench_Newton(n, c, abar, iter, pi_tol)
    % read default parameters
    mp0=trmodel.default_params;
    
    % Parameter adjustments (relative to defaults)
    mp0.ntypes = n;
    mp0.ncartypes=c;
    mp0.abar=abar;
    mp0.pnew = 100;
    [I, J] = ndgrid(0:n-1, 0:c-1);
    u_0 = 5 + 2 * (I + J) / (n + c);
    mp0.u_0 = u_0;
    
    % populate mp with remaining fields and update model parameter dependencies
    mp = trmodel.setparams(mp0);
    
    % model indexing
    s = trmodel.index(mp);
    
    % initialize prices where to solve bellman
    p = mp.pnew.*0.85.^(1:mp.abar-1)'; 
    
    % pre-compute age transition matrices
    F = trmodel.age_transition(mp, s);
    
    ap=dpsolver.setup;
    ap.pi_max = iter;
    ap.pi_tol = pi_tol;

    % No tolerance


    ev0 = zeros(s.ns,1); % initial guess on expected value function

    ev = zeros(s.ns, mp.ntypes);
    iter_n = zeros(1, mp.ntypes);


    for t  = 1:mp.ntypes
        % pre-compute utility
        util = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau
        
        % step 1: solve the DP problem for the price vector p for each consumer type
        bellman=@(ev) trmodel.bellman(mp, s, util, F, ev);
        
        [ev(:,t), P, dv, iter]= dpsolver.nk(bellman, ev0, ap);
        ev(:,t) = trmodel.bellman(mp, s, util, F, ev(:,t));
        iter_n(t) = iter.n;
    end
end


% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline
