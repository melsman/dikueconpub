% Test 4: demand and supply at initial prices p0 — compare to Futhark's diag_demand_supply_p0
%
% Futhark call:
%   echo "2i64 1i64 25i64 20i64 0f64 [0.1f64, 0.3f64]" | ./run_equilibrium diag_demand_supply_p0
%
% This script computes demand and supply separately at the same initial prices
% that Futhark uses (simple_prices with r=0.85, pnew=200), so we can detect
% any index shifts in deltaT.

addpath('matlabinclude');
addpath('autotrade');
addpath('..');

clear
clc

%% Setup — match Futhark parameters exactly
mp = setparams.default();
mp.ntypes   = 2;
mp.ncartypes = 1;
mp.lbl_cartypes = {' '};
mp.transcost = 0;
mp.mum = {0.1; 0.3};   % matches Futhark mum=[0.1, 0.3]
mp = trmodel.update_mp(mp);
s = trmodel.index(mp);

%% Initial prices — match Futhark's simple_prices(mp, r=0.85) with pnew=200
% Futhark: p[ct][a] = pnew * r^a for a=1..Ax-1 (used cars, 0-indexed age)
% Matlab:  p(s.ip{j}) = prices for car type j, ages 1..abar-1
p0 = 200 * 0.85 .^ (1:s.np)';   % s.np = abar_j0{1}-1 = 24

%% Initialize ev0 global (as done in equilibrium.solve_once)
global ev0;
if isempty(ev0) || numel(ev0) ~= mp.ntypes || numel(ev0{1}) ~= s.ns
    for t = 1:mp.ntypes
        ev0{t} = zeros(s.ns, 1);
    end
end

%% Evaluate edf at p0 to get solution components
[~, ~, sol] = equilibrium.edf(mp, s, p0);

%% Compute demand and supply per type — matching Futhark's demand_supply_tau
idx = [s.is.car_ex_clunker{:}];   % state indices for used cars (ages 1..abar-1)

demand_tau = zeros(s.np, mp.ntypes);
supply_tau = zeros(s.np, mp.ntypes);

for t = 1:mp.ntypes
    if mp.tw(t) > 0
        [~, ~, deltaT_t, ~] = trmodel.trade_transition(mp, s, sol.ccp_tau{t}, sol.ccp_scrap_tau{t});
        demand_tau(:, t) = deltaT_t(:, idx)' * sol.q_tau{t};
        supply_tau(:, t) = (1 - sol.ccp_tau{t}(idx, s.id.keep)) .* ...
                           (1 - sol.ccp_scrap_tau{t}(idx, :)) .* ...
                           sol.q_tau{t}(idx);
    end
end

%% Aggregate over types (weighted by mp.tw)
demand = demand_tau * mp.tw(:);
supply = supply_tau * mp.tw(:);

%% Print results for comparison with Futhark
fprintf('Initial prices p0 (ages 1..%d):\n', s.np);
disp(p0');

fprintf('\nDemand (ages 1..%d):\n', s.np);
disp(demand');

fprintf('\nSupply (ages 1..%d):\n', s.np);
disp(supply');

fprintf('\nExcess demand (demand - supply):\n');
disp((demand - supply)');

fprintf('\nDemand vs Supply table (age | demand | supply | ed):\n');
fprintf('%5s %12s %12s %12s\n', 'age', 'demand', 'supply', 'ed');
for a = 1:s.np
    fprintf('%5d %12.6f %12.6f %12.6f\n', a, demand(a), supply(a), demand(a)-supply(a));
end
