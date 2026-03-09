% This solves the bellman equation for a given set of prices and parameters

function [ev, ev_sa, ev_poly]=run_bellman_mul(n, c, abar, type, pri, futhark_u0)

    
    if nargin<4
        type = "all";
    end
    if nargin<5
        pri = true;
    end
    
    mp = setparams.default();

    mp.tw =  ones(n,1)/n;
    
    % Parameter adjustments (relative to defaults)
    mp.ncartypes=c;
    mp.abar_j0 ={abar};
    mp.ntypes=n;
    mp.pnew={100};
    mp.pnew_notax={100};
    mp.acc_0 = {-10};
    mp.u_0 = {6};
    mp.u_a = {-0.5};
    mp.u_a_sq = {0.0};
    mp.transcost = 0;
    mp.vat = 0;
    mp.cartax_lo =  0;   % registration tax (below kink, K_cartax_hi)
    mp.cartax_hi =  0;   % registration tax (above kink, K_cartax_hi)
    mp.tax_fuel  =  0;   % proportional fuel tax 
    mp.K_cartax_hi = 0;       % mp.K_cartax_hi before mp.cartax_hi tax kicks in
    mp.mum = {0.1};

    if nargin==6 & futhark_u0==true
        [I, J] = ndgrid(0:n-1, 0:c-1);
        u_0 = 5 + 2 * (I + J) / (n + c);
        mp.u_0 = num2cell(u_0);
    end

    
    % update model parameter dependencies
    mp = trmodel.update_mp(mp);

    if pri
        fprintf('%.6f\n', mp.u_0);
    end
    
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
    if pri
        ap.printfxp=1;  % show output
    end
    ev0 = zeros(s.ns,1); % initial guess on expected value function

    ev_sa   = zeros(s.ns, mp.ntypes);
    ev_poly = zeros(s.ns, mp.ntypes);

    ev = zeros(s.ns, mp.ntypes);

    for t  = 1:mp.ntypes
       if pri
        fprintf('Solving model for sonsumer type, tau=%d\n ', t);
       end
        
        % pre-compute utility
        util = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau
        
        % step 1: solve the DP problem for the price vector p for each consumer type
        bellman=@(ev) trmodel.bellman(mp, s, util, F, ev);
        
        if type~="poly"
            ev_sa(:,t)= dpsolver.sa(bellman, ev0, ap);
            ev(:,t) = ev_sa(:,t);
        end
        if type~="sa"
            ev_poly(:,t)= dpsolver.poly(bellman, ev0, ap, mp.bet);
            ev(:,t) = ev_poly(:,t);
        end
    end
end
