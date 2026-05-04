function [sol, mp, s] = bench_run_illustrations(n, c, abar, acc_0, transcost)
    mp = setparams.default(); % parameters used for illustration
    mp.ntypes = n;
    mp.tw = ones(n,1)/n;
    mp.ncartypes=c;
    mp.abar_j0 = {abar};
    mp.acc_0 = {acc_0};
    mp.transcost = transcost;
    mp = trmodel.update_mp(mp);
    s = trmodel.index(mp);

    sol=equilibrium.solve(mp, s);
end
