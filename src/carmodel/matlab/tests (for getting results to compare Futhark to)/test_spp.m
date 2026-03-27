% Test script for spp.bellman_spp, spp.solve_spp, and equilibrium.price_init

addpath('matlabinclude');
addpath('autotrade');
addpath('..');

clear
clc

%% Setup: same parameters as run_illustrations_test
mp = setparams.default();
mp.ntypes = 2;
mp.ncartypes = 1;
mp.lbl_cartypes = {' '};
mp.abar_j0 = {20, 20};
mp.transcost = 0;
mp.acc_0 = {-5};
mp = trmodel.update_mp(mp);

%% 1. Test bellman_spp: single Bellman iteration from zero
fprintf('=== bellman_spp test ===\n');
abar = spp.abar_max; % same abar used inside solve_spp
v0 = zeros(abar, 1);

% Test for type 1, car 1
[V1, P1, dV1, Eu1] = spp.bellman_spp(mp, v0, 1, 1);
fprintf('Type 1, Car 1:\n');
fprintf('  V (first 10): '); fprintf('%.6f ', V1(1:min(10,end))); fprintf('\n');
fprintf('  V (last 5):   '); fprintf('%.6f ', V1(end-4:end)); fprintf('\n');
fprintf('  P (first 10): '); fprintf('%d ', full(P1(1:min(10,end)))); fprintf('\n');
fprintf('  sum(P):        %d\n', full(sum(P1)));

% Test for type 2, car 1
[V2, P2, dV2, Eu2] = spp.bellman_spp(mp, v0, 2, 1);
fprintf('Type 2, Car 1:\n');
fprintf('  V (first 10): '); fprintf('%.6f ', V2(1:min(10,end))); fprintf('\n');
fprintf('  V (last 5):   '); fprintf('%.6f ', V2(end-4:end)); fprintf('\n');
fprintf('  P (first 10): '); fprintf('%d ', full(P2(1:min(10,end)))); fprintf('\n');
fprintf('  sum(P):        %d\n', full(sum(P2)));

%% 2. Test solve_spp: full solution
fprintf('\n=== solve_spp test ===\n');
[abar_spp, p_spp, w_spp] = spp.solve_spp(mp);

for t = 1:mp.ntypes
    for car = 1:mp.ncartypes
        fprintf('Type %d, Car %d:\n', t, car);
        fprintf('  abar_spp: %d\n', abar_spp{t, car});
        fprintf('  p_spp:    '); fprintf('%.6f ', p_spp{t, car}'); fprintf('\n');
        fprintf('  w_spp:    '); fprintf('%.6f ', w_spp{t, car}'); fprintf('\n');
    end
end

%% 3. Test equilibrium.price_init
fprintf('\n=== price_init test ===\n');
s = trmodel.index(mp);
p0 = equilibrium.price_init(mp, s);
fprintf('  p0 (length %d): ', length(p0));
fprintf('%.6f ', p0');
fprintf('\n');
