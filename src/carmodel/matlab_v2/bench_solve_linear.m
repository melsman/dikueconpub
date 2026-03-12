% This solves the bellman equation for a given set of prices and parameters


function [x, t] = bench_solve_linear(c, abar)
    f = @() bench_solve_linear_once(c, abar);
    t = timeit(f);
    x = bench_solve_linear_once(c, abar);
end

function x = bench_solve_linear_once(c, abar)
    n = 1;

    mp0 = trmodel.default_params;
    mp0.ntypes = n;
    mp0.ncartypes = c;
    mp0.abar = abar;
    mp0.pnew = 100;

    [I, J] = ndgrid(0:n-1, 0:c-1);
    u_0 = 5 + 2 * (I + J) / (n + c);
    mp0.u_0 = u_0;

    mp = trmodel.setparams(mp0);
    s = trmodel.index(mp);
    p = mp.pnew .* 0.85.^(1:mp.abar-1)';
    Ftr = trmodel.age_transition(mp, s);
    ev0 = zeros(s.ns,1);
    tau = 1;
    util = trmodel.utility(mp, s, tau, p);

    [ev1, ~, dV] = trmodel.bellman(mp, s, util, Ftr, ev0);
    A = speye(s.ns) - dV;
    x = ev0-A \ (ev0 - ev1);
end


% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline
