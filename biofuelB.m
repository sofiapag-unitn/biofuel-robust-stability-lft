clear;
clc;
close all;
%  Coupled Parametric LFT:  Biofuel Model B
%  Anchored at the MIDPOINT of the swept equilibrium box
fprintf('Biofuel Model B  (12 parameters: kinetic + equilibrium shift)\n');
fprintf('=============================================================\n\n');
% Nominal physical parameters
alpha_n = 0.66;
n_max   = 1.0;
delta_n = 0.91;
gamma_p = 0.14;
beta_p  = 0.66;
alpha_b = 0.10;
delta_b = 0.50;
alpha_p = 0.01;
k_p     = 0.20;
gamma_I = 60;
I       = 1000;
% repressor parameters (fixed -- not uncertain, so R_bar is constant)
alpha_R = 0.01;
beta_R  = 2.1;
k_R     = 10;
gamma_R = 1.8;
k_b     = 100;
rho = 0.20;
% Rational terms as anonymous functions (used by J0 and the finite differences)
phi_fun = @(gp, pb) gp ./ (pb + gp).^2;                         % phi_p
psi_fun = @(kp, kb, gR, Rb, bb) ...                             % psi_b = d(pdot)/d(b_i)
          kp .* Rb .* kb ./ ((Rb ./ (1 + kb.*bb) + gR).^2 .* (1 + kb.*bb).^2);
% Repressor steady state -- FIXED (alpha_R, beta_R, k_R, I, gamma_I are all
% fixed, so R_bar has no uncertain inputs left and never varies).
R_bar = (alpha_R + k_R*(I/(I + gamma_I))) / beta_R;
fprintf('Fixed repressor steady state:  R_bar = %.6f  (constant, not swept)\n\n', R_bar);
% Nominal equilibrium
model_ode = @(x) [
    x(1) * (alpha_n*(1 - x(1)/n_max) - delta_n*x(3) ...
            - alpha_n*x(2)/(x(2) + gamma_p));
    alpha_p + k_p * (1 / (R_bar/(1 + k_b*x(3)) + gamma_R)) - beta_p*x(2);
    alpha_b*x(1) - delta_b*x(2)*x(3)
];
opts  = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-13);
x_bar = fsolve(model_ode, [0.17; 0.17; 0.21], opts);
n_bar  = x_bar(1);     % nominal equilibrium (for reference)
p_bar  = x_bar(2);
bi_bar = x_bar(3);
fprintf('Nominal equilibrium:\n');
fprintf('  n_bar  = %.6f\n', n_bar);
fprintf('  p_bar  = %.6f\n', p_bar);
fprintf('  bi_bar = %.6f\n', bi_bar);
% Swept equilibrium box over the +/-20% parameter box (9 kinetic parameters
% only -- alpha_R, beta_R, k_R, alpha_p are fixed), then MIDPOINT anchor.
% extremes confirmed at the corners against a 50000-point interior check
n_bar_lo = 0.0724;  n_bar_hi = 0.2909;
p_bar_lo = 0.0957;  p_bar_hi = 0.2978;
b_bar_lo = 0.0751;  b_bar_hi = 0.4954;
n_ctr = 0.5*(n_bar_lo + n_bar_hi);   rho_n_bar = 0.5*(n_bar_hi - n_bar_lo)/n_ctr;
p_ctr = 0.5*(p_bar_lo + p_bar_hi);   rho_p_bar = 0.5*(p_bar_hi - p_bar_lo)/p_ctr;
b_ctr = 0.5*(b_bar_lo + b_bar_hi);   rho_b_bar = 0.5*(b_bar_hi - b_bar_lo)/b_ctr;
% Anchor equilibrium = midpoints (kinetic parameters stay at nominal)
n_bar  = n_ctr;
p_bar  = p_ctr;
bi_bar = b_ctr;
% R_bar is NOT re-anchored -- it stays at its fixed value from above.
fprintf('\nAnchoring at equilibrium-box midpoint:\n');
fprintf('  n_ctr  = %.6f   interval [%.4f, %.4f]\n', n_bar,  n_bar_lo, n_bar_hi);
fprintf('  p_ctr  = %.6f   interval [%.4f, %.4f]\n', p_bar,  p_bar_lo, p_bar_hi);
fprintf('  b_ctr  = %.6f   interval [%.4f, %.4f]\n', bi_bar, b_bar_lo, b_bar_hi);
fprintf('  R_bar  = %.6f   (fixed, no interval)\n', R_bar);
fprintf('Half-width radii (symmetric about midpoint):  rho_nbar=%.3f, rho_pbar=%.3f, rho_bbar=%.3f\n\n', ...
    rho_n_bar, rho_p_bar, rho_b_bar);
% Jacobian J0 at the midpoint anchor
phi_p = phi_fun(gamma_p, p_bar);
psi_b = psi_fun(k_p, k_b, gamma_R, R_bar, bi_bar);
J0 = zeros(3, 3);
J0(1,1) = -alpha_n * n_bar / n_max;
J0(1,2) = -alpha_n * n_bar * phi_p;
J0(1,3) = -delta_n * n_bar;
J0(2,2) = -beta_p;
J0(2,3) =  psi_b;
J0(3,1) =  alpha_b;
J0(3,2) = -delta_b * bi_bar;
J0(3,3) = -delta_b * p_bar;
fprintf('phi_p = %.4f,   psi_b = %.6f   (feedback entry J0(2,3))\n', phi_p, psi_b);
fprintf('Jacobian J0 at midpoint anchor:\n');
disp(J0);
eig_J0 = eig(J0);
fprintf('Eigenvalues of J0:\n');
disp(eig_J0);
fprintf('max Re eig(J0) = %+.4f   (stable if < 0)\n\n', max(real(eig_J0)));
%  per-parameter sensitivity matrices Jk
%  for nonlinear parameters Jk is computed via finite difference
tol_rank = 1e-9;
ek       = 1e-6;
fprintf('Sensitivity matrices Jk:\n\n');
% linear kinetic parameters: Jk(i,j) = rho * J0(i,j)
% alpha_n entries (1,1) and (1,2)
Jk_alpha_n      = zeros(3,3);
Jk_alpha_n(1,1) = rho * J0(1,1);
Jk_alpha_n(1,2) = rho * J0(1,2);
fprintf('  alpha_n  ->  Jk(1,1) = %+.5f,  Jk(1,2) = %+.5f  (linear, exact)\n', ...
    Jk_alpha_n(1,1), Jk_alpha_n(1,2));
% delta_n entry (1,3)
Jk_delta_n      = zeros(3,3);
Jk_delta_n(1,3) = rho * J0(1,3);
fprintf('  delta_n  ->  Jk(1,3) = %+.5f  (linear, exact)\n', Jk_delta_n(1,3));
% beta_p entry (2,2)
Jk_beta_p       = zeros(3,3);
Jk_beta_p(2,2)  = rho * J0(2,2);
fprintf('  beta_p   ->  Jk(2,2) = %+.5f  (linear, exact)\n', Jk_beta_p(2,2));
% alpha_b entry (3,1)
Jk_alpha_b      = zeros(3,3);
Jk_alpha_b(3,1) = rho * J0(3,1);
fprintf('  alpha_b  ->  Jk(3,1) = %+.5f  (linear, exact)\n', Jk_alpha_b(3,1));
% delta_b entries (3,2) and (3,3)
Jk_delta_b      = zeros(3,3);
Jk_delta_b(3,2) = rho * J0(3,2);
Jk_delta_b(3,3) = rho * J0(3,3);
fprintf('  delta_b  ->  Jk(3,2) = %+.5f,  Jk(3,3) = %+.5f  (linear, exact)\n', ...
    Jk_delta_b(3,2), Jk_delta_b(3,3));
% k_p entry (2,3)
Jk_k_p      = zeros(3,3);
Jk_k_p(2,3) = rho * J0(2,3);
fprintf('  k_p      ->  Jk(2,3) = %+.5f  (linear, exact)\n', Jk_k_p(2,3));
% nonlinear kinetic parameters: finite difference
% gamma_p entry (1,2)
gamma_p_p   = gamma_p * (1 + ek);
J_pert      = J0;
J_pert(1,2) = -alpha_n * n_bar * phi_fun(gamma_p_p, p_bar);
Jk_gamma_p  = rho * gamma_p * (J_pert - J0) / (gamma_p * ek);
fprintf('  gamma_p  ->  Jk(1,2) = %+.5f  (nonlinear, finite difference)\n', Jk_gamma_p(1,2));
% k_b entry (2,3)
k_b_p       = k_b * (1 + ek);
J_pert      = J0;
J_pert(2,3) = psi_fun(k_p, k_b_p, gamma_R, R_bar, bi_bar);
Jk_k_b      = rho * k_b * (J_pert - J0) / (k_b * ek);
fprintf('  k_b      ->  Jk(2,3) = %+.5f  (nonlinear, finite difference)\n', Jk_k_b(2,3));
% gamma_R entry (2,3)
gamma_R_p   = gamma_R * (1 + ek);
J_pert      = J0;
J_pert(2,3) = psi_fun(k_p, k_b, gamma_R_p, R_bar, bi_bar);
Jk_gamma_R  = rho * gamma_R * (J_pert - J0) / (gamma_R * ek);
fprintf('  gamma_R  ->  Jk(2,3) = %+.5f  (nonlinear, finite difference)\n', Jk_gamma_R(2,3));
% equilibrium-shift parameters (n_bar, p_bar, bi_bar only -- R_bar is fixed)
% n_bar: linear in all of row 1   Jk(1,j) = rho_n_bar * J0(1,j)
Jk_n_bar      = zeros(3,3);
Jk_n_bar(1,1) = rho_n_bar * J0(1,1);
Jk_n_bar(1,2) = rho_n_bar * J0(1,2);
Jk_n_bar(1,3) = rho_n_bar * J0(1,3);
fprintf('  n_bar    ->  Jk(1,1)=%+.5f, Jk(1,2)=%+.5f, Jk(1,3)=%+.5f  (linear, exact)\n', ...
    Jk_n_bar(1,1), Jk_n_bar(1,2), Jk_n_bar(1,3));
% p_bar: nonlinear in (1,2) linear in (3,3)
p_bar_p     = p_bar * (1 + ek);
J_pert      = J0;
J_pert(1,2) = -alpha_n * n_bar * phi_fun(gamma_p, p_bar_p);
J_pert(3,3) = -delta_b * p_bar_p;
Jk_p_bar    = rho_p_bar * p_bar * (J_pert - J0) / (p_bar * ek);
fprintf('  p_bar    ->  Jk(1,2) = %+.5f,  Jk(3,3) = %+.5f  (nonlinear, finite diff, rank 2)\n', ...
    Jk_p_bar(1,2), Jk_p_bar(3,3));
% b_bar: nonlinear in (2,3) linear in (3,2)
b_bar_p     = bi_bar * (1 + ek);
J_pert      = J0;
J_pert(2,3) = psi_fun(k_p, k_b, gamma_R, R_bar, b_bar_p);
J_pert(3,2) = -delta_b * b_bar_p;
Jk_b_bar    = rho_b_bar * bi_bar * (J_pert - J0) / (bi_bar * ek);
fprintf('  bi_bar   ->  Jk(2,3) = %+.5f,  Jk(3,2) = %+.5f  (nonlinear, finite diff, rank 2)\n\n', ...
    Jk_b_bar(2,3), Jk_b_bar(3,2));
%  Rank factorization  Jk = Lk * Rk   (SVD)
fprintf('Rank factorization Jk = Lk * Rk:\n');
[L1,  R1,  rk1 ] = rank1_factor(Jk_alpha_n, tol_rank);  fprintf('  alpha_n  rank = %d\n', rk1);
[L2,  R2,  rk2 ] = rank1_factor(Jk_delta_n, tol_rank);  fprintf('  delta_n  rank = %d\n', rk2);
[L3,  R3,  rk3 ] = rank1_factor(Jk_gamma_p, tol_rank);  fprintf('  gamma_p  rank = %d\n', rk3);
[L4,  R4,  rk4 ] = rank1_factor(Jk_beta_p,  tol_rank);  fprintf('  beta_p   rank = %d\n', rk4);
[L5,  R5,  rk5 ] = rank1_factor(Jk_alpha_b, tol_rank);  fprintf('  alpha_b  rank = %d\n', rk5);
[L6,  R6,  rk6 ] = rank1_factor(Jk_delta_b, tol_rank);  fprintf('  delta_b  rank = %d\n', rk6);
[L7,  R7,  rk7 ] = rank1_factor(Jk_k_p,     tol_rank);  fprintf('  k_p      rank = %d\n', rk7);
[L8,  R8,  rk8 ] = rank1_factor(Jk_k_b,     tol_rank);  fprintf('  k_b      rank = %d\n', rk8);
[L9,  R9,  rk9 ] = rank1_factor(Jk_gamma_R, tol_rank);  fprintf('  gamma_R  rank = %d\n', rk9);
[L10, R10, rk10] = rank1_factor(Jk_n_bar,   tol_rank);  fprintf('  n_bar    rank = %d\n', rk10);
[L11, R11, rk11] = rank1_factor(Jk_p_bar,   tol_rank);  fprintf('  p_bar    rank = %d  (repeated scalar)\n', rk11);
[L12, R12, rk12] = rank1_factor(Jk_b_bar,   tol_rank);  fprintf('  bi_bar   rank = %d  (repeated scalar)\n\n', rk12);
%  Stack L, R and build Delta block structure
%  Order: 9 kinetic, then n_bar, p_bar(x2), b_bar(x2).  (R_bar is fixed -- no channel.)
L_mat = [L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12];
R_mat = [R1; R2; R3; R4; R5; R6; R7; R8; R9; R10; R11; R12];
m = 12;   % physical uncertain parameters (9 kinetic + n_bar, p_bar, bi_bar)
param_ranks = [rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9, rk10, rk11, rk12];
q = sum(param_ranks);
fprintf('L matrix (3 x %d):\n', q);  disp(L_mat);
fprintf('R matrix (%d x 3):\n', q);  disp(R_mat);
%  Block structure
blk = zeros(m, 2);
for k = 1:m
    blk(k, 1) = -param_ranks(k);   % -1 for scalars, -2 for the two repeated blocks
    blk(k, 2) =  0;
end
display(blk);
%  Verification
Jk_sum = Jk_alpha_n + Jk_delta_n + Jk_gamma_p + Jk_beta_p + Jk_alpha_b + ...
         Jk_delta_b + Jk_k_p + Jk_k_b + Jk_gamma_R + Jk_n_bar + ...
         Jk_p_bar + Jk_b_bar;
max_error = max(abs(L_mat * eye(q) * R_mat - Jk_sum));
fprintf('Verification  max |L*I*R - sum(Jk)| = %.2e   (should be ~0)\n\n', max(max_error));
%  Pre-sweep conditioning check
n_states = 3;
epsilon  = 0.0005;
M_at_zero = R_mat * ((-J0) \ L_mat);
fprintf('--- Pre-sweep check ---\n');
fprintf('  norm(M) at w = 0          = %.4f\n', norm(M_at_zero));
fprintf('  Property 1 (mu)  at w = 0 = %.4f\n',  mu_upper_bound_property1(M_at_zero, blk));
fprintf('  Property 2 (gam) at w = 0 = %.4f\n',  mu_upper_bound_property2(M_at_zero, blk));
fprintf('  Perron           at w = 0 = %.4f\n',  gamma_perron_bound(M_at_zero, blk, epsilon));
fprintf('  Recursive        at w = 0 = %.4f\n',  mu_upper_bound_recursive(M_at_zero, blk));
[mb0, ~] = mussv(M_at_zero, blk);
fprintf('  mussv (repeated) at w = 0 = %.4f\n',  mb0(1));
fprintf('-----------------------\n\n');
%  Frequency sweep and bounds
w_grid = logspace(-3, 3, 15);
n_grid = length(w_grid);
mu_bound        = zeros(1, n_grid);
gamma_bound     = zeros(1, n_grid);
gamma_perron    = zeros(1, n_grid);
recursive_bound = zeros(1, n_grid);
mu_mussv        = zeros(1, n_grid);
fprintf('Computing bounds across %d frequencies...\n', n_grid);
for k = 1:n_grid
    w = w_grid(k);
    s = 1j * w;
    M = R_mat * ((s * eye(n_states) - J0) \ L_mat);
    mu_bound(k)        = mu_upper_bound_property1(M, blk);
    gamma_bound(k)     = mu_upper_bound_property2(M, blk);
    recursive_bound(k) = mu_upper_bound_recursive(M, blk);
    gamma_perron(k)    = gamma_perron_bound(M, blk, epsilon);
    [mu_b, ~]   = mussv(M, blk);
    mu_mussv(k) = mu_b(1);
end
fprintf('Done.\n\n');
%  Results
mu_max        = max(max(mu_bound),        mu_upper_bound_property1(M_at_zero, blk));
gamma_max     = max(max(gamma_bound),     mu_upper_bound_property2(M_at_zero, blk));
recursive_max = max(max(recursive_bound), mu_upper_bound_recursive(M_at_zero, blk));
perron_max    = max(max(gamma_perron),    gamma_perron_bound(M_at_zero, blk, epsilon));
mussv_max     = max(max(mu_mussv),        mb0(1));
fprintf('======== RESULTS  (Model B, 12 parameters) ========\n');
fprintf('Property 1 (mu)    sup = %.4f   ->  ||delta||_inf  < %.4f   \n', mu_max,        1/mu_max);
fprintf('Property 2 (gamma) sup = %.4f   ->  ||delta||_2    < %.4f   \n', gamma_max,     1/gamma_max);
fprintf('Recursive          sup = %.4f   ->  ||delta||_2    < %.4f   \n', recursive_max, 1/recursive_max);
fprintf('Perron             sup = %.4f   ->  ||delta||_2    < %.4f   \n', perron_max,  1/perron_max);
fprintf('mussv              sup = %.4f   ->  ||delta||_inf  < %.4f   \n', mussv_max, 1/mussv_max);
inf_radius      = 1 / mu_max;
two_norm_radius = 1 / gamma_max;
sqrt_m          = sqrt(m);
fprintf('\nBox coverage  (delta in [-1,1]^m,  m = %d physical parameters):\n', m);
if inf_radius >= 1
    fprintf('  inf-norm radius  = %.4f >= 1.0                  : box COVERED\n', inf_radius);
else
    fprintf('  inf-norm radius  = %.4f < 1.0                   : box NOT covered\n', inf_radius);
end
if two_norm_radius >= sqrt_m
    fprintf('  2-norm radius           = %.4f >= sqrt(%d) = %.4f   : box COVERED\n', two_norm_radius, m, sqrt_m);
else
    fprintf('  2-norm radius           = %.4f < sqrt(%d) = %.4f    : box NOT covered\n', two_norm_radius, m, sqrt_m);
end
%  Plots
figure('Color', 'w');
semilogx(w_grid, mu_bound, 'r--', 'LineWidth', 1.5);
hold on;
semilogx(w_grid, gamma_bound, 'b-', 'LineWidth', 1.5);
semilogx(w_grid, mu_mussv, 'k-', 'LineWidth', 1.5);
xlabel('\omega (rad/s)');  ylabel('Upper bounds');
title('Biofuel B (12 params): \mu vs \gamma upper bounds');
legend('\mu (Property 1)', '\gamma (Property 2)', '\mu (mussv, repeated)');
grid on;  box on;
figure('Color', 'w');
semilogx(w_grid, gamma_bound, 'r--', 'LineWidth', 1.5);
hold on;
semilogx(w_grid, gamma_perron, 'b-', 'LineWidth', 1.5);
semilogx(w_grid, recursive_bound, 'k-', 'LineWidth', 1.5);
xlabel('\omega (rad/s)');  ylabel('Upper bounds');
title('Biofuel B (12 params): \gamma bounds comparison');
legend('Property 2', 'Perron bound', 'Recursive bound');
grid on;  box on;
%  Local helper: SVD rank factorization  Jk = Lk * Rk
function [Lk, Rk, rk] = rank1_factor(Jk, tol_rank)
    [U, S, V] = svd(Jk);
    s  = diag(S);
    rk = sum(s > tol_rank * max(s(1), eps));
    Lk = U(:, 1:rk) * diag(sqrt(s(1:rk)));
    Rk = diag(sqrt(s(1:rk))) * V(:, 1:rk)';
end