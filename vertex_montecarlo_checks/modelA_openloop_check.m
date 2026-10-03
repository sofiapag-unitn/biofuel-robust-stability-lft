function modelA_openloop_check()
% Multiaffinity + grid check for Biofuel Model A
% sweeps the BDC diagonal channels D(1..7) 
clear;
clc;
% nominal parameter values
alpha_n = 0.66;    n_max   = 1.0;
gamma_p = 0.14;    delta_n = 0.91;
beta_p  = 0.66;    alpha_b = 0.10;
delta_b = 0.50;
alpha_p = 0.01;    k_p     = 0.20;
I       = 1000.0;  gamma_I = 60.0;

RHO = 0.20;        % +/-20% box

% The six uncertain parameters
%   1: alpha_n   2: delta_n   3: gamma_p   4: beta_p   5: alpha_b   6: delta_b
nomU  = [alpha_n, delta_n, gamma_p, beta_p, alpha_b, delta_b];
names = {'alpha_n','delta_n','gamma_p','beta_p','alpha_b','delta_b'};
coordNames = {'n_bar','p_bar','b_bar'};
Dnames = {'D1 an*n/nmax','D2 an*n*phi','D3 dn*n','D4 beta_p', ...
          'D5 alpha_b','D6 db*b','D7 db*p'};

% nominal check
d0 = [0, 0, 0, 0, 0, 0];
[n0, p0, b0, D0] = equilibrium_from_deltas(d0, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
fprintf('Nominal x_bar = [%.4f %.4f %.4f]  \n\n', n0, p0, b0);

% per-axis non-affinity (curvature / range; 0.00 = exactly affine)
fprintf('Per-axis non-affinity (curvature/range; 0.00 = exactly affine):\n');
fprintf('%9s | %8s %8s %8s\n', 'param', 'n_bar', 'p_bar', 'b_bar');
nLevels    = 11;
axisLevels = linspace(-1, 1, nLevels);
for k = 1:6
    nvals = zeros(nLevels, 1);
    pvals = zeros(nLevels, 1);
    bvals = zeros(nLevels, 1);
    for j = 1:nLevels
        d = [0, 0, 0, 0, 0, 0];
        d(k) = axisLevels(j);
        [nn, pp, bb] = equilibrium_from_deltas(d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
        nvals(j) = nn;
        pvals(j) = pp;
        bvals(j) = bb;
    end
    sn = curvature_score(nvals);
    sp = curvature_score(pvals);
    sb = curvature_score(bvals);
    fprintf('%9s | %8.3f %8.3f %8.3f\n', names{k}, sn, sp, sb);
end

% vertices (2^6) and a fine grid (7^6 points)
% scan_box sweeps every combination of the given gridValues across all parameters:
%   gridValues = [-1, 1]           the 64 vertices
%   gridValues = linspace(-1,1,7)  the grid
% It also returns the per-entry BDC D-box over the same samples
[vmin, vmax, Dvmin, Dvmax] = scan_box([-1, 1],            nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
[gmin, gmax, Dgmin, Dgmax] = scan_box(linspace(-1, 1, 7), nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);

% report (equilibrium coordinates)
fprintf('\nInterval over the +/-20%% box:\n');
fprintf('%6s | %8s %8s | %8s %8s | extremes at vertices?\n', ...
        'coord', 'vtx min', 'vtx max', 'grid min', 'grid max');
for c = 1:3
    atVertices = (gmin(c) >= vmin(c) - 1e-9) && (gmax(c) <= vmax(c) + 1e-9);
    if atVertices
        tag = 'YES';
    else
        tag = 'NO (interior extremum found)';
    end
    fprintf('%6s | %8.4f %8.4f | %8.4f %8.4f | %s\n', ...
            coordNames{c}, vmin(c), vmax(c), gmin(c), gmax(c), tag);
end

fprintf('\nCentering radii (centered at NOMINAL equilibrium):\n');
x0 = [n0, p0, b0];
for c = 1:3
    lo = gmin(c);  hi = gmax(c);
    r_abs = max(x0(c) - lo, hi - x0(c));   % worst-case one-sided deviation
    rho   = r_abs / x0(c);
    fprintf('  %-6s nominal=%.4f  rho=%.3f (%.1f%%)   [down %.1f%%, up %.1f%%]\n', ...
        coordNames{c}, x0(c), rho, 100*rho, ...
        100*(x0(c) - lo)/x0(c), 100*(hi - x0(c))/x0(c));
end

% BDC box: vertex vs interior test 
fprintf('\nBDC box over the +/-20%% box  (vertices 2^6=64  vs  7^6 grid):\n');
fprintf('%-13s | %9s | %9s %9s | %9s %9s | interior within vertices?\n', ...
        'entry', 'D_nom', 'vtx min', 'vtx max', 'grid min', 'grid max');
allYES = true;
for c = 1:7
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

% Final D-box = union of vertex and grid extremes 
Dmin = min(Dvmin, Dgmin);
Dmax = max(Dvmax, Dgmax);
fprintf('\nFinal D-box (union of vertex + grid) and percent about D_nom:\n');
fprintf('%-13s | %9s %9s | %8s %8s\n', 'entry', 'D_min', 'D_max', 'down%', 'up%');
for c = 1:7
    fprintf('%-13s | %9.4g %9.4g | %8.1f %8.1f\n', ...
            Dnames{c}, Dmin(c), Dmax(c), ...
            100*(1 - Dmin(c)/D0(c)), 100*(Dmax(c)/D0(c) - 1));
end
if allYES
    fprintf(['\n==> All D-entries take their extremes at the vertices \n']);
else
    fprintf(['\n==> An interior D-extremum was found \n']);
end

% Monte Carlo / probabilistic stability comparison 
% Builds the EXACT Jacobian J = B*diag(D)*C at sampled points and tests
% Hurwitz (max Re eig < 0)
%   (1) all 2^m vertices         
%   (1b) uniform interior sample 
%   (2) independent BDC D-box    
%   (3) radius sweep             -
B = zeros(3,7);
B(1,1)=-1; B(1,2)=-1; B(1,3)=-1;
B(2,4)=-1;
B(3,5)= 1; B(3,6)=-1; B(3,7)=-1;
C = zeros(7,3);
C(1,1)=1; C(2,2)=1; C(3,3)=1;
C(4,2)=1; C(5,1)=1; C(6,2)=1; C(7,3)=1;

m = numel(nomU);
fprintf('\n================ Monte Carlo stability comparison ================\n');
fprintf('(Reference from the mu analysis: BDC 1/mu = 1.5099 covered;  LFT 1/mu = 1.3599 covered)\n');

% (1) deterministic worst case over all 2^m vertices
[nUns_v, worstRe_v] = vertex_stability(nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C);
fprintf('\nVertices (all 2^%d = %d corners): %d unstable;  worst max Re eig = %+.4e\n', ...
        m, 2^m, nUns_v, worstRe_v);

% (1b) interior Monte Carlo (uniform in the box) 
Nmc = 50000;
[ns, Np, worstRe_p, nBad] = mc_physical(Nmc, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C);
fprintf(' MC (%d uniform draws): P(stable) = %.5f', Np, ns/Np);
if ns == Np
    fprintf('   [0 unstable : P(unstable) < %.2e @95%%]\n', 3/Np);
else
    fprintf('   [%d unstable]\n', Np - ns);
end
fprintf('   worst max Re eig = %+.4e', worstRe_p);
if nBad > 0, fprintf('   (%d non-physical draws skipped)', nBad); end
fprintf('\n');

% (2) BDC independent-D-box Monte Carlo
[nsb, Nb, worstRe_b] = mc_bdc_set(50000, B, C, Dmin, Dmax);
fprintf('BDC-set MC (%d independent D draws): P(stable) = %.5f   [%d unstable]; worst max Re eig = %+.4e\n', ...
        Nb, nsb/Nb, Nb - nsb, worstRe_b);

% (3) radius sweep: scale the box by alpha and estimate P(stable)
alphas = [0.50 0.80 1.00 1.20 1.3599 1.50 1.5099 1.70 2.00 2.50 3.00 3.50 4.00 4.50 4.8 4.9 5];
[al, Ps] = mc_radius_sweep(alphas, 5000, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C);
fprintf('\nRadius sweep (box scaled by alpha, 5000 draws each):\n');
fprintf('   %9s %12s\n', 'alpha', 'P(stable)');
for a = 1:numel(al)
    fprintf('   %9.4f %12.5f\n', al(a), Ps(a));
end

% (3b) SAME alphas, but VERTEX-based: check all 2^m corners of the scaled box
[alv, wRe, nUnsw, nParBad, nEqBad] = ...
    vertex_radius_sweep(alphas, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C);
fprintf('\nVERTEX-based radius sweep (worst over all 2^%d corners of the scaled box):\n', m);
fprintf('   %9s %16s %10s %9s %8s\n', 'alpha', 'worst maxRe', '#unstable', '#par<=0', '#eq<=0');
for a = 1:numel(alv)
    fprintf('   %9.4f %16.4e %10d %9d %8d\n', alv(a), wRe(a), nUnsw(a), nParBad(a), nEqBad(a));
end
fprintf('==================================================================\n');

end


% =====================================================================
function [xmin, xmax, Dmin, Dmax] = scan_box(gridValues, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max)
% Evaluate the equilibrium and the BDC D-vector at every combination of the
% per-axis values and return the min/max of each coordinate and each D entry.
%   gridValues = [-1, 1]           the 128 vertices
%   gridValues = linspace(-1,1,7)  the grid

xmin = [ Inf,  Inf,  Inf];
xmax = [-Inf, -Inf, -Inf];
Dmin =  inf(7, 1);
Dmax = -inf(7, 1);
for v1 = gridValues
  for v2 = gridValues
    for v3 = gridValues
      for v4 = gridValues
        for v5 = gridValues
          for v6 = gridValues
                d = [v1, v2, v3, v4, v5, v6];
                [nn, pp, bb, D] = equilibrium_from_deltas(d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
                x = [nn, pp, bb];
                for c = 1:3
                   if x(c) < xmin(c)
                      xmin(c) = x(c);
                   end
                   if x(c) > xmax(c)
                     xmax(c) = x(c);
                   end
                end
                Dmin = min(Dmin, D);
                Dmax = max(Dmax, D);
          end
        end
      end
    end
  end
end

end


% =====================================================================
function [n, p, b, D] = equilibrium_from_deltas(d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max)
% Build the perturbed parameters from the delta vector d (length 7, in [-1,1]),
% evaluate the closed-form open-loop equilibrium, and assemble the BDC
% diagonal channels D(1..7) 

alpha_n = nomU(1)*(1 + RHO*d(1));
delta_n = nomU(2)*(1 + RHO*d(2));
gamma_p = nomU(3)*(1 + RHO*d(3));
beta_p  = nomU(4)*(1 + RHO*d(4));
alpha_b = nomU(5)*(1 + RHO*d(5));
delta_b = nomU(6)*(1 + RHO*d(6));


% Pump equation is decoupled
C = alpha_p + k_p*I/(I + gamma_I);
p = C/beta_p;

% Cell density
numerator   = alpha_n*gamma_p*n_max*delta_b*p;
denominator = (p + gamma_p)*(alpha_n*delta_b*p + delta_n*alpha_b*n_max);
n = numerator/denominator;

% biofuel
b = alpha_b*n/(delta_b*p);

% BDC diagonal channels 
phi_p = gamma_p/(p + gamma_p)^2;
D = [ alpha_n*n/n_max;      % 1
      alpha_n*n*phi_p;      % 2
      delta_n*n;            % 3
      beta_p;               % 4 (pure parameter)
      alpha_b;              % 5 (pure parameter)
      delta_b*b;            % 6
      delta_b*p ];          % 7

end

% =====================================================================
function score = curvature_score(y)
% Largest discrete second difference relative to the total range.
% 0 means the values lie on a straight line
rng = max(y) - min(y);
if rng > 1e-12
    secondDiff = diff(y, 2);
    score = max(abs(secondDiff))/rng;
else
    score = 0.0;
end

end

% =====================================================================
function [nUns, worstRe] = vertex_stability(nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C)
% Exact Jacobian J = B*diag(D)*C at all 2^m parameter corners; count unstable.
m = numel(nomU);
nUns = 0; worstRe = -inf;
for idx = 0:(2^m - 1)
    bits       = bitget(idx, 1:m);     % 0/1 across the m parameters
    d          = 2*bits - 1;           % -1 / +1 corner
    [~,~,~, D] = equilibrium_from_deltas(d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
    mr = max(real(eig(B*diag(D)*C)));
    if mr >= 0, nUns = nUns + 1; end
    if mr > worstRe, worstRe = mr; end
end
end

% =====================================================================
function [nStable, N, worstRe, nBad] = mc_physical(N, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C)
%  uniform draws in the +/-20% box, exact Jacobian, Hurwitz test.
rng(1);
m = numel(nomU);
nStable = 0; nBad = 0; worstRe = -inf;
for i = 1:N
    d = 2*rand(1, m) - 1;
    [n, p, b, D] = equilibrium_from_deltas(d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
    if n <= 0 || p <= 0 || b <= 0 || any(~isfinite(D))
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
function [alphas, Pstab] = mc_radius_sweep(alphas, N, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C)
% Scale the box by alpha (param = nom*(1 + alpha*RHO*d)); estimate P(stable).
% Non-physical equilibria (negative states) are skipped, not counted.
rng(3);
m = numel(nomU);
Pstab = zeros(size(alphas));
for a = 1:numel(alphas)
    al = alphas(a); ns = 0; nv = 0;
    for i = 1:N
        d = 2*rand(1, m) - 1;
        [n, p, b, D] = equilibrium_from_deltas(al*d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
        if n <= 0 || p <= 0 || b <= 0 || any(~isfinite(D)), continue; end
        nv = nv + 1;
        if max(real(eig(B*diag(D)*C))) < 0, ns = ns + 1; end
    end
    if nv > 0, Pstab(a) = ns/nv; else, Pstab(a) = NaN; end
end
end

% =====================================================================
function [alphas, worstRe, nUns, nParBad, nEqBad] = vertex_radius_sweep(alphas, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max, B, C)
% At each alpha, scan all 2^m corners of the scaled box and report the worst
% max Re eig over corners that are PHYSICALLY VALID (positive parameters AND a
% positive equilibrium).  

%   nParBad - corners where a perturbed parameter went <= 0
%   nEqBad  - corners with positive parameters but a non-positive equilibrium

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
            nParBad(a) = nParBad(a) + 1;
            continue;
        end
        [n, p, b, D] = equilibrium_from_deltas(al*d, nomU, RHO,alpha_p, k_p, I, gamma_I, n_max);
        if n <= 0 || p <= 0 || b <= 0 || any(~isfinite(D))
            nEqBad(a) = nEqBad(a) + 1;
            continue;
        end
        mr = max(real(eig(B*diag(D)*C)));
        if mr > worstRe(a), worstRe(a) = mr; end
        if mr >= 0, nUns(a) = nUns(a) + 1; end
    end
end
end
