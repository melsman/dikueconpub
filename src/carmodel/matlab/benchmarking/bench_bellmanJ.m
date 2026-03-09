% This solves the bellman equation for a given set of prices and parameters


function [ev, ccp, dev] = bench_bellmanJ(n, c, abar, N)
    % read default parameters
    mp0=setparams.default();
    
    % Parameter adjustments (relative to defaults)
    mp0.ncartypes=c;
    mp0.abar_j0 ={abar};
    mp0.ntypes=n;
    mp0.pnew={100};
    mp0.pnew_notax={100};
    mp0.acc_0 = {-10};
    mp0.u_0 = {6};
    mp0.u_a = {-0.5};
    mp0.u_a_sq = {0.0};
    mp0.transcost = 0;
    mp0.vat = 0;
    mp0.cartax_lo =  0;   % registration tax (below kink, K_cartax_hi)
    mp0.cartax_hi =  0;   % registration tax (above kink, K_cartax_hi)
    mp0.tax_fuel  =  0;   % proportional fuel tax 
    mp0.K_cartax_hi = 0;       % mp.K_cartax_hi before mp.cartax_hi tax kicks in
    mp0.mum = {0.1};
    
    % populate mp with remaining fields and update model parameter dependencies
    mp = trmodel.setparams(mp0);
    
    % model indexing
    s = trmodel.index(mp);
    
    % initialize prices where to solve bellman
    function [sp] = simple_prices(mp, r,s)
        sp=nan(s.np,1);
        for j=1:mp.ncartypes;
            i=1:(mp.abar_j0{j}-1);
            sp(s.ip{j}) = mp.pnew_notax{j}.*r.^i;
        end
    end
    p = simple_prices(mp,0.85,s); 
    
    % pre-compute age transition matrices
    F = trmodel.age_transition(mp, s);
    
    ap=dpsolver.setup;
    ap.printfxp=1;  % show output
    
    t=1;
    
    % pre-compute utility
    util = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau
    
    % step 1: solve the DP problem for the price vector p for each consumer type
    bellman=@(ev) trmodel.bellman(mp, s, util, F, ev);
    
    ev0=zeros(s.ns, 1); % initial guess on expected value function
    for i=1:N
        [ev0(:,1), ccp, dev] = bellman (ev0(:,1));
    end
    
    ev = ev0(:,1);
end


% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline
