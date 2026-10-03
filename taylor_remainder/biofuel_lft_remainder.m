clear;
clc;
close all;
%  Coupled Parametric LFT:  Biofuel Models A and B
%  SYMBOLIC TOTAL-DERIVATIVE construction, PLUS a Taylor-form remainder

%    WAY A: residual sampling
%      Re-solve the TRUE equilibrium and TRUE Jacobian at box corners and
%      interior Monte-Carlo draws and take the worst-case entrywise deviation from the affine surrogate

%    WAY B: analytic curvature via the recursive total-derivative operator
%      The SAME total-derivative operator applied once more to tot_k gives the
%      Hessian term that Taylor's theorem needs for the remainder:

% SELECT MODEL HERE 
model = 'A';        % 'A'  ->  Model A only  (states [n;p;b],   q = 13)
% 'B'  ->  Model B only  (states [n;R;p;b], q = 27)
% 'both' -> run A then B

rho     = 0.20;
Nmc     = 300;         % interior Monte-Carlo draws for WAY A
rng(1);

runA = any(strcmpi(model, {'A','both'}));
runB = any(strcmpi(model, {'B','both'}));

% all symbols (declared once so either block can run alone)
syms n p b R real
syms alpha_n delta_n gamma_p beta_p alpha_b delta_b real      % shared kinetic
syms alpha_p k_p gamma_I I n_max real                         % fixed in A (k_p uncertain in B)
syms k_b gamma_R alpha_R k_R beta_R real                      % Model B extras

%  MODEL A
if runA
    fprintf('################  MODEL A  ################\n\n');
    xA     = [n; p; b];
    thetaA = [alpha_n; delta_n; gamma_p; beta_p; alpha_b; delta_b];
    fA = [ n*( alpha_n*(1 - n/n_max) - delta_n*b - alpha_n*p/(p + gamma_p) );
           alpha_p + k_p*(I/(I + gamma_I)) - beta_p*p;
           alpha_b*n - delta_b*p*b ];
    parsA = [alpha_n; delta_n; gamma_p; beta_p; alpha_b; delta_b; alpha_p; k_p; gamma_I; I; n_max];
    valsA = [0.66;    0.91;    0.14;    0.66;   0.10;    0.50;    0.01;    0.20; 60;     1000; 1.0];
% closed-form equilibrium 
    pA = (alpha_p + k_p*(I/(I+gamma_I)))/beta_p;
    nA = (alpha_n*gamma_p/(pA+gamma_p)) / (alpha_n/n_max + delta_n*alpha_b/(delta_b*pA));
    bA = alpha_b*nA/(delta_b*pA);
    xbarA = [nA; pA; bA];                                  % SYMBOLIC equilibrium

    [L_A, R_A, blk_A, J0_A, qA, ranksA, ctxA] = build_lft_sym(fA, xA, thetaA, xbarA, parsA, valsA, rho);

    fprintf('\n--- Model A remainder (WAY A: residual sampling) ---\n');
    RmatA_A = remainder_sampling(ctxA, Nmc);
    fprintf('||R_A||_fro = %.4e , ||R_A||_2 = %.4e\n', norm(RmatA_A,'fro'), norm(RmatA_A,2));

    fprintf('\n--- Model A remainder (WAY B: analytic Hessian) ---\n');
    [HblkA, RmatA_B] = remainder_hessian(ctxA);
    fprintf('||R_B||_fro = %.4e , ||R_B||_2 = %.4e\n', norm(RmatA_B,'fro'), norm(RmatA_B,2));

    RmatA = max(RmatA_A, RmatA_B);        % combined, conservative envelope
    fprintf('\nCombined remainder bound |R| (entrywise max of A,B):\n'); disp(RmatA);

    [L_A2, R_A2, blk_A2] = augment_with_remainder(L_A, R_A, blk_A, RmatA);
    qA2 = size(L_A2,2);

    fprintf('\n=== Model A : baseline (first-order local) ===\n');
    run_mu(L_A, R_A, blk_A, J0_A, qA, ranksA, 'A (baseline)');
    fprintf('\n=== Model A : with remainder (global box) ===\n');
    run_mu(L_A2, R_A2, blk_A2, J0_A, qA2, ranksA, 'A (with remainder)');
end

%  MODEL B
if runB
    fprintf('\n\n################  MODEL B   ################\n\n');
    xB     = [n; R; p; b];
    thetaB = [alpha_n; delta_n; gamma_p; beta_p; alpha_b; delta_b; k_p; k_b; gamma_R];
    Reff   = R/(1 + k_b*b);
    fB = [ n*( alpha_n*(1 - n/n_max) - delta_n*b - alpha_n*p/(p + gamma_p) );
           alpha_R + k_R*(I/(I + gamma_I)) - beta_R*R;
           alpha_p + k_p/(Reff + gamma_R) - beta_p*p;
           alpha_b*n - delta_b*p*b ];
    parsB = [alpha_n; delta_n; gamma_p; beta_p; alpha_b; delta_b; k_p; k_b; gamma_R; ...
             alpha_p; gamma_I; I; n_max; alpha_R; k_R; beta_R];
    valsB = [0.66; 0.91; 0.14; 0.66; 0.10; 0.50; 0.20; 100; 1.8; ...
             0.01; 60; 1000; 1.0; 0.01; 10; 2.1];
% equilibrium: R* closed form; n*(p),b*(p) reduce fB(3) to a scalar g(p)=0
    RB     = (alpha_R + k_R*(I/(I+gamma_I)))/beta_R;
    n_of_p = (alpha_n*gamma_p/(p+gamma_p)) / (alpha_n/n_max + delta_n*alpha_b/(delta_b*p));
    b_of_p = subs(alpha_b*n/(delta_b*p), n, n_of_p);
    gpB    = subs(fB(3), [R, b], [RB, b_of_p]);
    gp_nom = subs(gpB, parsB, valsB);
    [gnum,~] = numden(gp_nom);  cp = sym2poly(gnum);
    rts = roots(cp);  rts = rts(abs(imag(rts))<1e-9 & real(rts)>0);
    psB = real(min(rts));
    RsB = double(subs(RB, parsB, valsB));
    nsB = double(subs(subs(n_of_p, parsB, valsB), p, psB));
    bsB = double(subs(subs(b_of_p, parsB, valsB), p, psB));
    xbarB_num = [nsB; RsB; psB; bsB];                      % NUMERIC equilibrium

    [L_B, R_B, blk_B, J0_B, qB, ranksB, ctxB] = build_lft_sym(fB, xB, thetaB, xbarB_num, parsB, valsB, rho);

    fprintf('\n--- Model B remainder (WAY A: residual sampling) ---\n');
    RmatB_A = remainder_sampling(ctxB, Nmc);
    fprintf('||R_A||_fro = %.4e , ||R_A||_2 = %.4e\n', norm(RmatB_A,'fro'), norm(RmatB_A,2));

    fprintf('\n--- Model B remainder (WAY B: analytic Hessian) ---\n');
    [HblkB, RmatB_B] = remainder_hessian(ctxB);
    fprintf('||R_B||_fro = %.4e , ||R_B||_2 = %.4e\n', norm(RmatB_B,'fro'), norm(RmatB_B,2));

    RmatB = max(RmatB_A, RmatB_B);
    fprintf('\nCombined remainder bound |R| (entrywise max of A,B):\n'); disp(RmatB);

    [L_B2, R_B2, blk_B2] = augment_with_remainder(L_B, R_B, blk_B, RmatB);
    qB2 = size(L_B2,2);

    fprintf('\n=== Model B : baseline (first-order local) ===\n');
    run_mu(L_B, R_B, blk_B, J0_B, qB, ranksB, 'B (baseline)');
    fprintf('\n=== Model B : with remainder (global box) ===\n');
    run_mu(L_B2, R_B2, blk_B2, J0_B, qB2, ranksB, 'B (with remainder)');
end


%  build_lft_sym : symbolic total-derivative LFT, evaluated at xbar
%  Now also returns Wblk (per-channel full W_k, before SVD splitting) and
%  a ctx struct bundling everything WAY A / WAY B need 
function [L_mat, R_mat, blk, J0, q, ranks, ctx] = build_lft_sym(f, x, theta, xbar, pars, vals, rho)
    nx = numel(x);  m = numel(theta);
    Fx    = jacobian(f, x);                 % symbolic state Jacobian
    Fth   = jacobian(f, theta);             % symbolic parameter Jacobian
    Sens  = -Fx \ Fth;                      % symbolic dxbar/dtheta  (IFT)
    dFx_dx = cell(1,nx);
for i = 1:nx, dFx_dx{i} = diff(Fx, x(i)); end
if isnumeric(xbar), xb = xbar(:); else, xb = double(subs(xbar, pars, vals)); end
    subVars = [x; pars];  subVals = [xb; vals];
    J0 = double(subs(Fx, subVars, subVals));
    eJ = eig(J0);
    fprintf('Equilibrium  xbar = [%s]\n', strtrim(sprintf('%.5f ', xb)));
    fprintf('J0 =\n'); disp(J0);
    fprintf('eig(J0): %s   (max Re = %+.4f)\n\n', ...
            strtrim(sprintf('%.4f%+.4fi  ', [real(eJ) imag(eJ)].')), max(real(eJ)));
    theta_nom = double(subs(theta, pars, vals));
    tol = 1e-9;  Lblk = cell(1,m); Rblk = cell(1,m); Wblk = cell(1,m);
    totblk = cell(1,m); ranks = zeros(m,1); W_sum = zeros(nx);
    fprintf('Per-channel  W_k = rho*thetabar_k*(dFx/dtheta_k)_total  ->  W_k = Lk*Rk:\n');
for k = 1:m
        tot = diff(Fx, theta(k));                          % explicit
for i = 1:nx, tot = tot + dFx_dx{i}*Sens(i,k); end % + implicit (symbolic)
        totblk{k} = tot;                                    % stash for WAY B
        Wk = rho * theta_nom(k) * double(subs(tot, subVars, subVals));
        Wblk{k} = Wk;                                        % stash for WAY A
        W_sum = W_sum + Wk;
        [U,S,V] = svd(Wk);  s = diag(S);  rk = sum(s > tol*max(s(1),eps));
        Lblk{k} = U(:,1:rk)*diag(sqrt(s(1:rk)));
        Rblk{k} = diag(sqrt(s(1:rk)))*V(:,1:rk)';
        ranks(k) = rk;
        fprintf('  %-8s rank = %d\n', char(theta(k)), rk);
end
    L_mat = [Lblk{:}];  R_mat = cat(1, Rblk{:});  q = sum(ranks);
    blk = [-ranks, zeros(m,1)];                            % repeated REAL scalars
    fprintf('\nL (%dx%d), R (%dx%d), Delta = blkdiag(delta_k I_rk), q = %d\n', ...
            size(L_mat,1),size(L_mat,2),size(R_mat,1),size(R_mat,2), q);
    fprintf('Verification  max |L*I*R - sum W_k| = %.2e\n', ...
            max(abs(L_mat*eye(q)*R_mat - W_sum), [], 'all'));

    % bundle everything WAY A / WAY B will need
    ctx.f = f;  ctx.x = x;  ctx.theta = theta;  ctx.pars = pars;  ctx.vals = vals;
    ctx.Fx = Fx;  ctx.Sens = Sens;  ctx.dFx_dx = dFx_dx;  ctx.tot = totblk;
    ctx.Wblk = Wblk;  ctx.xb = xb;  ctx.theta_nom = theta_nom;  ctx.rho = rho;
end


%  remainder_sampling : WAY A 
%  Compares the TRUE nonlinear Jacobian (equilibrium re-solved by fsolve,
%  starting from the nominal xbar) against the affine surrogate
%  J0 + sum_k delta_k*Wk, over box corners plus interior Monte-Carlo draws,
%  and returns the worst-case entrywise |residual|.
function Rmat = remainder_sampling(ctx, Nmc)
    m  = numel(ctx.theta);
    nx = numel(ctx.xb);
    theta_nom = ctx.theta_nom(:)';  rho = ctx.rho;

    % sanity check: theta must be the first m entries of pars, in order
if ~isequal(ctx.theta(:), ctx.pars(1:m))
        error('remainder_sampling assumes theta == pars(1:m); adjust pars ordering.');
end

    fnum  = matlabFunction(ctx.f,  'Vars', {ctx.x, ctx.pars});
    Fxnum = matlabFunction(ctx.Fx, 'Vars', {ctx.x, ctx.pars});
    J0    = double(subs(ctx.Fx, [ctx.x; ctx.pars], [ctx.xb; ctx.vals]));

    % corner set (mirrors the empirical-radius script's convention)
if m <= 13
        Sgn = 2*(dec2bin(0:2^m-1) - '0') - 1;
else
        Sgn = sign(randn(8192, m));  Sgn(Sgn==0) = 1;
end
    Dsamples = [Sgn; 2*rand(Nmc,m)-1];        % corners + interior MC draws

    opts = optimoptions('fsolve','Display','off');
    Rmat = zeros(nx,nx);
    nkept = 0;  nskip = 0;
for s = 1:size(Dsamples,1)
        delta = Dsamples(s,:);
        pars_s = ctx.vals(:);
        pars_s(1:m) = theta_nom .* (1 + rho*delta);
if any(pars_s(1:m) <= 0), nskip = nskip+1; continue; end   % non-physical parameters

        [xbar_s, ~, flag] = fsolve(@(xx) fnum(xx, pars_s), ctx.xb, opts);
if flag <= 0 || any(xbar_s <= 0) || any(~isfinite(xbar_s))
            nskip = nskip+1; continue;                       % no valid equilibrium
end

        Jtrue   = Fxnum(xbar_s, pars_s);
        Jaffine = J0;
for k = 1:m
            Jaffine = Jaffine + delta(k)*ctx.Wblk{k};
end
        Rmat = max(Rmat, abs(Jtrue - Jaffine));
        nkept = nkept + 1;
end
    fprintf('  sampled %d points (%d kept, %d skipped non-physical/non-converged)\n', ...
            size(Dsamples,1), nkept, nskip);
end


%  remainder_hessian : WAY B -- analytic curvature via the recursive
%  total-derivative operator.
%     H_jk = D_thetaj[ tot_k ]
%          = d(tot_k)/dtheta_j  +  sum_i d(tot_k)/dx_i * Sens(i,j)
%  evaluated at the nominal equilibrium, scaled to delta-space, then
%  folded into a conservative entrywise bound using
%     delta_k^2       in [0,1]   (diagonal terms)
%     delta_j*delta_k in [-1,1]  (off-diagonal terms, j~=k)
function [Hblk, Rmat] = remainder_hessian(ctx)
    m  = numel(ctx.theta);
    nx = numel(ctx.xb);
    x = ctx.x;  theta = ctx.theta;  pars = ctx.pars;
    subVars = [x; pars];  subVals = [ctx.xb; ctx.vals];
    theta_nom = ctx.theta_nom;  rho = ctx.rho;

    Hblk = cell(m,m);
for k = 1:m
        totk = ctx.tot{k};                                  % symbolic dFx/dtheta_k (total)
for j = 1:m
            Hjk = diff(totk, theta(j));                     % explicit part
for i = 1:nx
                Hjk = Hjk + diff(totk, x(i))*ctx.Sens(i,j);  % implicit, chain rule through xbar
end
            Hnum = double(subs(Hjk, subVars, subVals));
            Hblk{j,k} = rho^2 * theta_nom(j) * theta_nom(k) * Hnum;   % delta-space scaling
end
end

    % symmetry check (H_jk should equal H_kj up to numerical noise)
    asym = 0;
for j = 1:m
for k = 1:m
            asym = max(asym, max(abs(Hblk{j,k}(:) - Hblk{k,j}(:))));
end
end
    fprintf('  symmetry check  max |H_jk - H_kj| = %.2e\n', asym);

    Rmat = zeros(nx,nx);
for k = 1:m
        Rmat = Rmat + 0.5*abs(Hblk{k,k});                   % delta_k^2 in [0,1]
end
for j = 1:m-1
for k = j+1:m
            Rmat = Rmat + abs(Hblk{j,k});                    % delta_j*delta_k in [-1,1]
end
end
end


%  augment_with_remainder : fold an entrywise remainder bound Rmat into the
%  existing (L,R,blk) LFT.

function [L2, R2, blk2] = augment_with_remainder(L_mat, R_mat, blk, Rmat)
    nx  = size(Rmat,1);
    tol = 1e-10*max(abs(Rmat(:)));
    Lext = {};  Rext = {};  blkext = zeros(0,2);
for i = 1:nx
for j = 1:nx
if abs(Rmat(i,j)) > tol
                ei = zeros(nx,1); ei(i) = 1;
                ej = zeros(1,nx); ej(j) = 1;
                Lext{end+1} = ei;                       % nx x 1
                Rext{end+1} = Rmat(i,j) * ej;            % 1  x nx
                blkext = [blkext; -1, 0];                 % repeated real scalar, size 1
end
end
end
    L2 = [L_mat, Lext{:}];
    R2 = [R_mat; cat(1, Rext{:})];
    blk2 = [blk; blkext];
    fprintf('Remainder folded in as %d independent rank-1 real-scalar blocks (one per nonzero entry of Rmat)\n', numel(Lext));
    fprintf('q: %d -> %d ,  total blocks: %d -> %d\n', size(L_mat,2), size(L2,2), size(blk,1), size(blk2,1));
end


%  run_mu : build M(s) = R*(sI - J0)^{-1}*L and sweep mu / bounds
function run_mu(L_mat, R_mat, blk, J0, q, ranks, tag)
    nx = size(J0,1);  m = size(blk,1);   % number of blocks 
    have = @(f) exist(f,'file')==2 || exist(f,'builtin')==5;
    w = logspace(-3,3,5);  ng = numel(w);
    mu = nan(1,ng); p1 = nan(1,ng); p2 = nan(1,ng); rc = nan(1,ng); pr = nan(1,ng);
    eps0 = 5e-4;
for j = 1:ng
        M = R_mat * ((1j*w(j)*eye(nx) - J0) \ L_mat);
if have('mussv'),                   mb = mussv(M, blk);                 mu(j) = mb(1); end
if have('mu_upper_bound_property1'), p1(j) = mu_upper_bound_property1(M, blk);          end
if have('mu_upper_bound_property2'), p2(j) = mu_upper_bound_property2(M, blk);          end
if have('mu_upper_bound_recursive'), rc(j) = mu_upper_bound_recursive(M, blk);          end
if have('gamma_perron_bound'),       pr(j) = gamma_perron_bound(M, blk, eps0);          end
end
    fprintf('\n--- Model %s results (M is %dx%d, %d channels) ---\n', tag, q, q, m);
if any(~isnan(mu)), fprintf('mussv     sup mu = %.4f -> ||delta||_inf < %.4f\n', max(mu), 1/max(mu)); end
if any(~isnan(p1)), fprintf('Property1 sup    = %.4f -> ||delta||_inf < %.4f\n', max(p1), 1/max(p1)); end
if any(~isnan(p2)), fprintf('Property2 sup    = %.4f -> ||delta||_2   < %.4f\n', max(p2), 1/max(p2)); end
if any(~isnan(rc)), fprintf('Recursive sup    = %.4f -> ||delta||_2   < %.4f\n', max(rc), 1/max(rc)); end
if any(~isnan(pr)), fprintf('Perron    sup    = %.4f -> ||delta||_2   < %.4f\n', max(pr), 1/max(pr)); end
if any(~isnan(mu))
        figure('Color','w'); semilogx(w, mu, 'k-', 'LineWidth', 1.6); hold on;
if any(~isnan(p1)), semilogx(w, p1, 'r--'); end
if any(~isnan(p2)), semilogx(w, p2, 'b-.'); end
        xlabel('\omega (rad/s)'); ylabel('\mu upper bound');
        title(sprintf('Biofuel %s: \\mu vs \\omega', tag));
        grid on; box on;
end
end
