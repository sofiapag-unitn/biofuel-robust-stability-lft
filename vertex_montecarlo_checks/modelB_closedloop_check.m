function modelB_closedloop_check()
% Multiaffinity + grid/sample check for Biofuel Model B
% also sweeps the BDC diagonal channels D(1..8)
clear;
clc;
%  The equilibrium depends on 9 uncertain physical parameters so
% the dense grid of the Model A check is replaced by a large RANDOM INTERIOR SAMPLE
% nominal parameter values
alpha_n = 0.66;    n_max   = 1.0;
gamma_p = 0.14;    delta_n = 0.91;
beta_p  = 0.66;    alpha_b = 0.10;
delta_b = 0.50;
alpha_p = 0.01;    k_p     = 0.20;
I       = 1000.0;  gamma_I = 60.0;
alpha_R = 0.01;    beta_R  = 2.1;
k_R     = 10.0;    gamma_R = 1.8;
k_b     = 100.0;
RHO = 0.20;        % +/-20% box
%   1: alpha_n  2: delta_n  3: gamma_p  4: beta_p  5: alpha_b  6: delta_b
%   7: k_p      8: k_b      9: gamma_R
% Fixed (not uncertain): alpha_R, beta_R, k_R, alpha_p -> R_bar is therefore constant.
nomU  = [alpha_n, delta_n, gamma_p, beta_p, alpha_b, delta_b, ...
         k_p, k_b, gamma_R];
names = {'alpha_n','delta_n','gamma_p','beta_p','alpha_b','delta_b', ...
'k_p','k_b','gamma_R'};
coordNames = {'n_bar','p_bar','b_bar','R_bar'};
Dnames = {'D1 an*n/nmax','D2 an*n*phi','D3 dn*n','D4 beta_p', ...
'D5 psi_b','D6 alpha_b','D7 db*b','D8 db*p'};
m      = numel(nomU);
nCoord = 4;
qD     = 8;
% fsolve settings
opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12);
x0   = [0.17; 0.17; 0.21];   % initial guess near nominal equilibrium
% nominal check
d0 = zeros(1, m);
[xb, D0] = equilibrium_from_deltas(d0, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
fprintf(['Nominal x_bar = [%.4f %.4f %.4f %.4f]  \n'], ...
         xb(1), xb(2), xb(3), xb(4));
% per-axis non-affinity
fprintf('Per-axis non-affinity (curvature/range; 0.00 = exactly affine):\n');
fprintf('%9s | %8s %8s %8s %8s\n', 'param', coordNames{:});
nLevels    = 11;
axisLevels = linspace(-1, 1, nLevels);
for k = 1:m
    vals = zeros(nLevels, nCoord);
for j = 1:nLevels
        d    = zeros(1, m);
        d(k) = axisLevels(j);
        vals(j,:) = equilibrium_from_deltas(d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
end
    s = zeros(1, nCoord);
for c = 1:nCoord
        s(c) = curvature_score(vals(:,c));
end
    fprintf('%9s | %8.3f %8.3f %8.3f %8.3f\n', names{k}, s(1), s(2), s(3), s(4));
end
% vertices (2^9 = 512 corners) and a random interior sample
% Both scans now also return the per-entry BDC D-box over the same samples.
Nint = 50000;
[vmin, vmax, Dvmin, Dvmax] = scan_vertices(nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, m, nCoord, qD);
[gmin, gmax, Dgmin, Dgmax] = scan_interior(Nint, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, m, nCoord, qD);
% report (equilibrium coordinates)
fprintf('\nInterval over the +/-20%% box  (vertices: 2^%d = %d corners; interior: %d random):\n', ...
        m, 2^m, Nint);
fprintf('%6s | %8s %8s | %8s %8s | extremes at vertices?\n', ...
'coord', 'vtx min', 'vtx max', 'rand min', 'rand max');
for c = 1:nCoord
    atVertices = (gmin(c) >= vmin(c) - 1e-9) && (gmax(c) <= vmax(c) + 1e-9);
if atVertices
        tag = 'YES';
else
        tag = 'NO (interior extremum found)';
end
    fprintf('%6s | %8.4f %8.4f | %8.4f %8.4f | %s\n', ...
            coordNames{c}, vmin(c), vmax(c), gmin(c), gmax(c), tag);
end
% centering radii
fprintf('\nCentering radii (from the box interval):\n');
fprintf('%6s | %10s %8s | %10s %8s\n', ...
'coord', 'nom-center', 'rho_nom', 'midpoint', 'rho_mid');
for c = 1:nCoord
    lo  = min(vmin(c), gmin(c));
    hi  = max(vmax(c), gmax(c));
    nom = xb(c);
    rho_nom = max(nom - lo, hi - nom) / nom;
    mid     = 0.5*(hi + lo);
    rho_mid = 0.5*(hi - lo) / mid;
    fprintf('%6s | %10.4f %8.3f | %10.4f %8.3f\n', ...
            coordNames{c}, nom, rho_nom, mid, rho_mid);
end
% BDC D-box: vertex vs interior test
fprintf('\nBDC D-box over the +/-20%% box  (vertices 2^%d=%d  vs  %d random interior):\n', ...
        m, 2^m, Nint);
fprintf('%-13s | %9s | %9s %9s | %9s %9s | interior within vertices?\n', ...
'entry', 'D_nom', 'vtx min', 'vtx max', 'rand min', 'rand max');
allYES = true;
for c = 1:qD
    within = (Dgmin(c) >= Dvmin(c) - 1e-9) && (Dgmax(c) <= Dvmax(c) + 1e-9);
if within
        tag = 'YES';
else
        tag = 'NO (interior extremum found)';
        allYES = false;
end
    fprintf('%-13s | %9.4g | %9.4g %9.4g | %9.4g %9.4g | %s\n', ...
            Dnames{c}, D0(c), Dvmin(c), Dvmax(c), Dgmin(c), Dgmax(c), tag);
end
% Final D-box = union of vertex and interior extremes
Dmin = min(Dvmin, Dgmin);
Dmax = max(Dvmax, Dgmax);
fprintf('\nFinal D-box (union of vertex + interior) and percent about D_nom:\n');
fprintf('%-13s | %9s %9s | %8s %8s\n', 'entry', 'D_min', 'D_max', 'down%', 'up%');
for c = 1:qD
    fprintf('%-13s | %9.4g %9.4g | %8.1f %8.1f\n', ...
            Dnames{c}, Dmin(c), Dmax(c), ...
            100*(1 - Dmin(c)/D0(c)), 100*(Dmax(c)/D0(c) - 1));
end
% Monte Carlo / probabilistic stability comparison
% Builds the EXACT Jacobian J = B*diag(D)*C at sampled points and tests
% Hurwitz (max Re eig < 0)
%   (1) all 2^m vertices
%   (1b) uniform physical sample
%   (2) independent BDC D-box
%   (3) radius sweep
B = zeros(3,8);
B(1,1)=-1; B(1,2)=-1; B(1,3)=-1;
B(2,4)=-1; B(2,5)= 1;
B(3,6)= 1; B(3,7)=-1; B(3,8)=-1;
C = zeros(8,3);
C(1,1)=1; C(2,2)=1; C(3,3)=1;
C(4,2)=1; C(5,3)=1;
C(6,1)=1; C(7,2)=1; C(8,3)=1;
fprintf('\n================ Monte Carlo stability comparison ================\n');
fprintf('(Reference from the mu analysis: BDC 1/mu = <UPDATE_AFTER_BDC_RERUN> ...;  LFT 1/mu = <UPDATE> covered)\n');
% (1) deterministic worst case over all 2^m vertices
[nUns_v, worstRe_v] = vertex_stability(nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, B, C);
fprintf('\nVertices (all 2^%d = %d corners): %d unstable;  worst max Re eig = %+.4e\n', ...
        m, 2^m, nUns_v, worstRe_v);
% (1b) physical Monte Carlo (uniform in the box)
Nmc = 20000;
[ns, Np, worstRe_p, nBad] = mc_physical(Nmc, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, B, C);
fprintf('Physical MC (%d uniform draws): P(stable) = %.5f', Np, ns/Np);
fprintf('   worst max Re eig = %+.4e', worstRe_p);
if nBad > 0, fprintf('   (%d non-physical draws skipped)', nBad); end
fprintf('\n');
% (2) BDC independent-D-box Monte Carlo (no equilibrium solve needed)
[nsb, Nb, worstRe_b] = mc_bdc_set(50000, B, C, Dmin, Dmax);
fprintf('BDC-set MC (%d independent D draws): P(stable) = %.5f   [%d unstable]; worst max Re eig = %+.4e\n', ...
        Nb, nsb/Nb, Nb - nsb, worstRe_b);
% (3) radius sweep: scale the box by alpha and estimate P(stable)
alphas = [0.50 0.80 1.00 1.10 1.20 1.30 1.40 1.50 1.75 2.00];
[al, Ps] = mc_radius_sweep(alphas, 2000, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, B, C);
fprintf('\nRadius sweep (box scaled by alpha, 2000 draws each):\n');
fprintf('   %9s %12s\n', 'alpha', 'P(stable)');
for a = 1:numel(al)
    fprintf('   %9.4f %12.5f\n', al(a), Ps(a));
end
% (3b) SAME alphas, but VERTEX-based: check all 2^m corners of the scaled box.
guesses = {[0.17;0.17;0.21], [0.05;0.05;0.05], [0.30;0.40;0.80], [0.10;0.10;0.10]};
[alv, wRe, nUnsw, nParBad, nEqBad] = ...
    vertex_radius_sweep(alphas, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, guesses, B, C);
fprintf('\nVERTEX-based radius sweep (worst over all 2^%d corners of the scaled box):\n', m);
fprintf('   %9s %16s %10s %9s %8s\n', 'alpha', 'worst maxRe', '#unstable', '#par<=0', '#eq<=0');
for a = 1:numel(alv)
    fprintf('   %9.4f %16.4e %10d %9d %8d\n', alv(a), wRe(a), nUnsw(a), nParBad(a), nEqBad(a));
end
fprintf('==================================================================\n');
end
% =====================================================================
function [xmin, xmax, Dmin, Dmax] = scan_vertices(nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, m, nCoord, qD)
% Evaluate the equilibrium and the BDC D-vector at all 2^9 = 512 corners.
gridValues = [-1, 1];      % each parameter sits at its -20% or +20% corner
xmin =  Inf(1, nCoord);
xmax = -Inf(1, nCoord);
Dmin =  inf(qD, 1);
Dmax = -inf(qD, 1);
for v1 = gridValues
for v2 = gridValues
for v3 = gridValues
for v4 = gridValues
for v5 = gridValues
for v6 = gridValues
for v7 = gridValues
for v8 = gridValues
for v9 = gridValues
              d = [v1, v2, v3, v4, v5, v6, v7, v8, v9];
              [x, D] = equilibrium_from_deltas(d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
              xmin = min(xmin, x);
              xmax = max(xmax, x);
              Dmin = min(Dmin, D);
              Dmax = max(Dmax, D);
end
end
end
end
end
end
end
end
end
end
% =====================================================================
function [xmin, xmax, Dmin, Dmax] = scan_interior(N, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, m, nCoord, qD)
% Draw N uniform random points from the interior of [-1,1]^m and return the
% min/max of each coordinate and each BDC D entry
rng(0);                       % fix the seed so the run is reproducible
xmin =  Inf(1, nCoord);
xmax = -Inf(1, nCoord);
Dmin =  inf(qD, 1);
Dmax = -inf(qD, 1);
for i = 1:N
    d = 2*rand(1, m) - 1;     % one random point inside the box: 9 values in [-1,1]
    [x, D] = equilibrium_from_deltas(d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
    xmin = min(xmin, x);
    xmax = max(xmax, x);
    Dmin = min(Dmin, D);
    Dmax = max(Dmax, D);
end
end
% =====================================================================
function [x, D] = equilibrium_from_deltas(d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0)
% Build the perturbed parameters from the delta vector d, solve the equilibrium,
% and assemble the BDC diagonal channels D(1..8).
% alpha_p, alpha_R, beta_R, k_R are fixed (not perturbed) -> R_bar is constant.
alpha_n = nomU(1)*(1 + RHO*d(1));
delta_n = nomU(2)*(1 + RHO*d(2));
gamma_p = nomU(3)*(1 + RHO*d(3));
beta_p  = nomU(4)*(1 + RHO*d(4));
alpha_b = nomU(5)*(1 + RHO*d(5));
delta_b = nomU(6)*(1 + RHO*d(6));
k_p     = nomU(7)*(1 + RHO*d(7));
k_b     = nomU(8)*(1 + RHO*d(8));
gamma_R = nomU(9)*(1 + RHO*d(9));
% Repressor decouples (constant now: no perturbed inputs)
R_bar = (alpha_R + k_R*I/(I + gamma_I)) / beta_R;
% Coupled 3-state steady state
F = @(z) [
    z(1)*(alpha_n*(1 - z(1)/n_max) - delta_n*z(3) - alpha_n*z(2)/(z(2) + gamma_p));
    alpha_p + k_p*(1/(R_bar/(1 + k_b*z(3)) + gamma_R)) - beta_p*z(2);
    alpha_b*z(1) - delta_b*z(2)*z(3)
];
z = fsolve(F, x0, opts);
x = [z(1), z(2), z(3), R_bar];
% BDC diagonal channels
phi_p   = gamma_p/(z(2) + gamma_p)^2;
denom_R = R_bar/(1 + k_b*z(3)) + gamma_R;
psi     = (k_p*R_bar*k_b)/(denom_R^2 * (1 + k_b*z(3))^2);
D = [ alpha_n*z(1)/n_max;    % 1
      alpha_n*z(1)*phi_p;    % 2
      delta_n*z(1);          % 3
      beta_p;                % 4 (pure parameter)
      psi;                   % 5  = J0(2,3) repressor feedback
      alpha_b;               % 6 (pure parameter)
      delta_b*z(3);          % 7
      delta_b*z(2) ];        % 8
end
% =====================================================================
function score = curvature_score(y)
% Largest discrete second difference relative to the total range
% 0 means the values lie on a straight line.
rng_y = max(y) - min(y);
if rng_y > 1e-12
    secondDiff = diff(y, 2);
    score = max(abs(secondDiff))/rng_y;
else
    score = 0.0;
end
end
% =====================================================================
function [nUns, worstRe] = vertex_stability(nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, B, C)
% Exact Jacobian J = B*diag(D)*C at all 2^m parameter corners; count unstable.
m = numel(nomU);
nUns = 0; worstRe = -inf;
for idx = 0:(2^m - 1)
    bits   = bitget(idx, 1:m);     % 0/1 across the m parameters
    d      = 2*bits - 1;           % -1 / +1 corner
    [~, D] = equilibrium_from_deltas(d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
    mr = max(real(eig(B*diag(D)*C)));
if mr >= 0, nUns = nUns + 1; end
if mr > worstRe, worstRe = mr; end
end
end
% =====================================================================
function [nStable, N, worstRe, nBad] = mc_physical(N, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, B, C)
% Ground-truth MC: uniform draws in the +/-20% box, exact Jacobian, Hurwitz test.
rng(1);
m = numel(nomU);
nStable = 0; nBad = 0; worstRe = -inf;
for i = 1:N
    d = 2*rand(1, m) - 1;
    [x, D] = equilibrium_from_deltas(d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
if any(x(1:3) <= 0) || any(~isfinite(D))
        nBad = nBad + 1;
continue;
end
    mr = max(real(eig(B*diag(D)*C)));
if mr < 0, nStable = nStable + 1; end
if mr > worstRe, worstRe = mr; end
end
end
% =====================================================================
function [nStable, N, worstRe] = mc_bdc_set(N, B, C, Dmin, Dmax)
% MC over the INDEPENDENT BDC D-box: each entry uniform in [Dmin,Dmax].
% No equilibrium solve
rng(2);
q = numel(Dmin);
nStable = 0; worstRe = -inf;
for i = 1:N
    D  = Dmin + (Dmax - Dmin).*rand(q, 1);
    mr = max(real(eig(B*diag(D)*C)));
if mr < 0, nStable = nStable + 1; end
if mr > worstRe, worstRe = mr; end
end
end
% =====================================================================
function [alphas, Pstab] = mc_radius_sweep(alphas, N, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0, B, C)
% Scale the box by alpha (param = nom*(1 + alpha*RHO*d)); estimate P(stable).
% Non-physical equilibria (negative states) are skipped, not counted.
rng(3);
m = numel(nomU);
Pstab = zeros(size(alphas));
for a = 1:numel(alphas)
    al = alphas(a); ns = 0; nv = 0;
for i = 1:N
        d = 2*rand(1, m) - 1;
        [x, D] = equilibrium_from_deltas(al*d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, x0);
if any(x(1:3) <= 0) || any(~isfinite(D)), continue; end
        nv = nv + 1;
if max(real(eig(B*diag(D)*C))) < 0, ns = ns + 1; end
end
if nv > 0, Pstab(a) = ns/nv; else, Pstab(a) = NaN; end
end
end
% =====================================================================
function [alphas, worstRe, nUns, nParBad, nEqBad] = vertex_radius_sweep(alphas, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, guesses, B, C)
% At each alpha, scan all 2^m corners of the scaled box and report the worst
% max Re eig over corners that are PHYSICALLY VALID (positive parameters
%   nParBad - corners where a perturbed parameter went <= 0
%   nEqBad  - corners with positive parameters but no positive equilibrium
% The equilibrium uses a multi-start solve (first positive root over several
% guesses), so a corner is only called non-physical if NO guess yields one.
m = numel(nomU);
worstRe = -inf(size(alphas));
nUns    = zeros(size(alphas));
nParBad = zeros(size(alphas));
nEqBad  = zeros(size(alphas));
for a = 1:numel(alphas)
    al = alphas(a);
for idx = 0:(2^m - 1)
        bits = bitget(idx, 1:m);
        d    = 2*bits - 1;                       % +/-1 corner
        pars = nomU .* (1 + RHO*al*d);           % perturbed parameters
if any(pars <= 0)
            nParBad(a) = nParBad(a) + 1;         % parameters must stay positive
continue;
end
        xok = false; D = [];
for g = 1:numel(guesses)                 % multi-start equilibrium
            [x, Dg] = equilibrium_from_deltas(al*d, nomU, RHO, alpha_p, alpha_R, beta_R, k_R, n_max, I, gamma_I, opts, guesses{g});
if all(x(1:3) > 0) && all(isfinite(Dg))
                D = Dg; xok = true; break;
end
end
if ~xok
            nEqBad(a) = nEqBad(a) + 1;           % no positive equilibrium found
continue;
end
        mr = max(real(eig(B*diag(D)*C)));
if mr > worstRe(a), worstRe(a) = mr; end
if mr >= 0, nUns(a) = nUns(a) + 1; end
end
end
end
