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

ccp_tau=cell(mp.ntypes,1);
ctp_tau=cell(mp.ntypes,1);
q_tau=cell(mp.ntypes,1);
delta_tau=cell(mp.ntypes,1); 
deltaK_tau=cell(mp.ntypes,1); 
deltaT_tau=cell(mp.ntypes,1); 
delta_scrap=cell(mp.ntypes,1); 
util_tau=cell(mp.ntypes,1); 
ev_scrap_tau=cell(mp.ntypes,1); 
ccp_scrap_tau=cell(mp.ntypes,1); 

% pre-compute age transition matrices
F = trmodel.age_transition(mp, s);

ap=dpsolver.setup;
ap.printfxp=2;  % show output

nm = 1

t = 0;
daccprob = [];

for t=1:nm;
	fprintf('tau=%d\n ', t)

	% pre-compute utility
	[util_tau{t}, ev_scrap_tau{t}, ccp_scrap_tau{t}] = trmodel.utility(mp, s, t, p);  % update the current period, post-trade payoffs stored in the global util_tau

    % step 1: solve the DP problem for the price vector p for each consumer type
    bellman=@(ev) trmodel.bellman(mp, s, util_tau{t}, F, ev);
    ev_t = ev0(:, t);
    ev_poly(:,t)= dpsolver.poly(bellman, ev_t, ap, mp.bet);

    [ev_t, ccp_t] = bellman (ev_poly(:,t));
    ccp_tau{t} = ccp_t;  % store the choice probabilities for each type 

    % trade, keep and scrap transition matrices
    [delta_tau{t}, deltaK_tau{t}, deltaT_tau{t}, delta_scrap{t}]=trmodel.trade_transition(mp, s, ccp_tau{t}, ccp_scrap_tau{t});  

    % state transition matrix
    ctp_tau{t}=delta_tau{t}*F.notrade; 
    
    [du, dccp_scrap] = g_trmodel.du_dp(mp, s, 1, ccp_scrap_tau{t});
    dbellman = g_trmodel.dbellman(mp, s, ev_t, ccp_t, ccp_scrap_tau{t}, ctp_tau{t}, delta_tau{t}, du, daccprob, 'prices');

    dev=(eye(s.ns)-mp.bet*ctp_tau{t})\dbellman;

    dv = g_trmodel.dv(mp, s, F, dev, du, daccprob, ev_t, 'prices');
    
    % calculate the gradient of the conditional choice probabilities
    dlogccp=g_trmodel.dlogccp_dp(mp, s, ccp_t, dv);

    % calculate the gradient of the conditional choice probabilities
    dccp=g_trmodel.dccp_dp(mp, s, ccp_t, dlogccp); 

    dctp=g_trmodel.dctp_dp(mp, s, F.notrade, dccp, delta_tau{t}, daccprob);
    


    % calculate the type-specific car distributions
    [q_tau{t}, dq] =ergodic(ctp_tau{t},dctp); 

    %ded=g_trmodel.ded(mp, s, ccp_tau, dccp, q_tau{t}, dq, ccp_scrap_tau{t}, dccp_scrap);

end

%disp(dbellman)
%disp(dev);
%disp(dctp);
%disp(du);
% disp(ctp_tau{t});
% disp(dctp);
disp(dq);
disp(ccp_t);
disp(dccp);
disp(dctp);
disp(dccp_scrap);

dpsell=zeros(s.ns,s.np); 
dpsell([s.is.car_ex_clunker{:}],:)=eye(s.np)*(1-mp.ptc_sale);
disp(dpsell);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline
