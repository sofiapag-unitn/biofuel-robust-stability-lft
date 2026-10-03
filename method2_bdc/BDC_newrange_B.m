% Robust stability analysis via mu and gamma upper bounds
% using the BDC Decomposition
% Evaluated over a per-entry asymmetric D-box swept from the +/-20% parameter box
% Midpoint-anchored (Delta in [-1,1])

clear;
clc;
close all;

% Select the model to analyze
%model_selection = 'Biofuel_A';
model_selection = 'Biofuel_B';

fprintf('Running Analysis for: %s \n\n', model_selection);

switch model_selection
    
    case 'Biofuel_A'
        alpha_n = 0.66; alpha_p = 0.01; alpha_b = 0.10;
        beta_p  = 0.66; delta_n = 0.91; delta_b = 0.50;
        gamma_p = 0.14; gamma_I = 60;   k_p     = 0.20;
        n_max   = 1;    I       = 1000;
        
        p_bar = (alpha_p + k_p * I / (I + gamma_I)) / beta_p;
        num = alpha_n * gamma_p / (p_bar + gamma_p);
        den = alpha_n / n_max + delta_n * alpha_b / (delta_b * p_bar);
        n_bar = num / den;
        bi_bar = alpha_b * n_bar / (delta_b * p_bar);
        
        phi_p = gamma_p / (p_bar + gamma_p)^2;
        
        q = 7;          
        n = 3;   
        
        B = zeros(n, q);
        B(1, 1) = -1;  B(1, 2) = -1;  B(1, 3) = -1;  
        B(2, 4) = -1;  
        B(3, 5) =  1;  B(3, 6) = -1;  B(3, 7) = -1;  
        
        C = zeros(q, n);
        C(1, 1) = 1;  C(2, 2) = 1;  C(3, 3) = 1;  
        C(4, 2) = 1;  C(5, 1) = 1;  C(6, 2) = 1;  C(7, 3) = 1;  
        
        D_nom = zeros(q, 1);
        D_nom(1) = alpha_n * n_bar / n_max;   
        D_nom(2) = alpha_n * n_bar * phi_p;   
        D_nom(3) = delta_n * n_bar;           
        D_nom(4) = beta_p;                    
        D_nom(5) = alpha_b;                   
        D_nom(6) = delta_b * bi_bar;          
        D_nom(7) = delta_b * p_bar;           

    case 'Biofuel_B'
        alpha_n = 0.66; alpha_p = 0.01; alpha_b = 0.10; alpha_R = 0.01;
        beta_R  = 2.1;  beta_p  = 0.66; 
        delta_n = 0.91; delta_b = 0.50;
        gamma_p = 0.14; gamma_I = 60;   gamma_R = 1.8;
        k_R     = 10;   k_p     = 0.20; k_b     = 100;
        n_max   = 1;    I       = 1000;
        
        R_bar = (alpha_R + k_R * (I / (I + gamma_I))) / beta_R;
        
        model_ode = @(x) [
            x(1) * (alpha_n * (1 - x(1)/n_max) - delta_n * x(3) - (alpha_n * x(2))/(x(2) + gamma_p));
            alpha_p + k_p * (1 / (R_bar / (1 + k_b * x(3)) + gamma_R)) - beta_p * x(2);
            alpha_b * x(1) - delta_b * x(2) * x(3)
        ];
        opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-13);
        x_bar = fsolve(model_ode, [0.17; 0.17; 0.21], opts);
        n_bar = x_bar(1); p_bar = x_bar(2); bi_bar = x_bar(3);
        
        phi_p = gamma_p / (p_bar + gamma_p)^2;
        denom_R = R_bar / (1 + k_b * bi_bar) + gamma_R;
        psi_bi = (k_p * R_bar * k_b) / (denom_R^2 * (1 + k_b * bi_bar)^2);
        
        q = 8;
        n = 3;
        
        B = zeros(n, q);
        B(1, 1) = -1;  B(1, 2) = -1;  B(1, 3) = -1; 
        B(2, 4) = -1;  B(2, 5) =  1; 
        B(3, 6) =  1;  B(3, 7) = -1;  B(3, 8) = -1; 
        
        C = zeros(q, n);
        C(1, 1) = 1;  C(2, 2) = 1;  C(3, 3) = 1; 
        C(4, 2) = 1;  C(5, 3) = 1; 
        C(6, 1) = 1;  C(7, 2) = 1;  C(8, 3) = 1; 
        
        D_nom = zeros(q, 1);
        D_nom(1) = alpha_n * n_bar / n_max;
        D_nom(2) = alpha_n * n_bar * phi_p;
        D_nom(3) = delta_n * n_bar;
        D_nom(4) = beta_p;
        D_nom(5) = psi_bi;
        D_nom(6) = alpha_b;
        D_nom(7) = delta_b * bi_bar;
        D_nom(8) = delta_b * p_bar;

    otherwise
        error('Invalid model selection.');
end

% Robustness evaluation

% Define the parameter box bounds
rho = 0.20;   % 20% around nominal 

% Diagonal scaling 
d_sqrt = sqrt(D_nom);

B_scaled     = B * diag(d_sqrt);
C_scaled     = diag(d_sqrt) * C;
D_nom_scaled = ones(q, 1);

% Per-entry D box, swept over the +/-rho PARAMETER box so that the shift of the
% equilibrium is accounted for (each D entry is evaluated directly at every sample)
switch model_selection
case 'Biofuel_A'
   D_lo = [0.04572; 0.02148; 0.07954; 0.528; 0.08; 0.02323; 0.1003];
   D_hi = [0.2134; 0.2044; 0.2528; 0.792; 0.12; 0.1107; 0.2258];
case 'Biofuel_B'
   D_lo = [0.03824; 0.03556; 0.07467; 0.528; 0.004814; 0.08; 0.03476; 0.03856];
   D_hi = [0.2304;  0.521;   0.2606;  0.792; 0.3677;    0.12; 0.263;   0.1746];
end

if ~isempty(D_lo)
    D_min_scaled = D_lo ./ D_nom;
    D_max_scaled = D_hi ./ D_nom;

    % per-entry rho readout
    fprintf('Per-entry D box (swept over the +/-%.0f%% parameter box):\n', 100*rho);
    fprintf(' entry |    D_nom    |  down%%     up%%\n');
    for i = 1:q
        fprintf('  D(%d)  | %10.4g | %7.1f  %7.1f\n', ...
            i, D_nom(i), 100*(1 - D_min_scaled(i)), 100*(D_max_scaled(i) - 1));
    end
    fprintf('\n');
else
    % fallback: symmetric +/- rho on D
    D_min_scaled = (1 - rho) * D_nom_scaled;
    D_max_scaled = (1 + rho) * D_nom_scaled;
end

% Anchor the BDC at the box MIDPOINT (centered, Delta in [-1,1]) 
D_ctr = 0.5*(D_min_scaled + D_max_scaled);   % box center 
W     = 0.5*(D_max_scaled - D_min_scaled);   % per-entry half-widths

J_nom = B_scaled * diag(D_nom_scaled) * C_scaled;   % nominal Jacobian (sanity check)
J_min = B_scaled * diag(D_min_scaled) * C_scaled;   % lower-corner Jacobian (reference)
J0    = B_scaled * diag(D_ctr)        * C_scaled;   % anchor Jacobian = midpoint

fprintf('max Re eig at D_min = %+.3e\n', max(real(eig(J_min))));
fprintf('max Re eig at D_nom = %+.3e\n', max(real(eig(J_nom))));
fprintf('max Re eig at D_ctr = %+.3e   (anchor; Hurwitz if < 0)\n', max(real(eig(J0))));

fprintf('Jacobian J0 (at D_ctr = box midpoint):\n');
disp(J0);
fprintf('Eigenvalues of J0:\n');
disp(eig(J0));

% Build M(s) directly around the midpoint anchor:  M(s) = C (sI - J0)^{-1} B diag(W)
M_func = @(s) C_scaled * ((s*eye(n) - J0) \ (B_scaled * diag(W)));


% Frequency sweep and bounds
w_grid = logspace(-3, 3, 15);
n_grid = length(w_grid);

mu_bound        = zeros(1, n_grid);
gamma_bound     = zeros(1, n_grid);
gamma_perron    = zeros(1, n_grid);
recursive_bound = zeros(1, n_grid);
mu_mussv        = zeros(1, n_grid);

% Block structure for mussv (q real repeated-scalar blocks of size 1)
blk = zeros(q, 2);
for i = 1:q
    blk(i, 1) = -1;
    blk(i, 2) = 0;
end

% conditioning check
M1  = M_func(1j*w_grid(1));
eJ0 = eig(J0);

fprintf('\n--- Pre-sweep check (w = %.1e) ---\n', w_grid(1));
fprintf('  norm(M)        = %.3e\n', norm(M1));
fprintf('  cond(J0)       = %.3e\n', cond(J0));     
fprintf('  max Re eig(J0) = %+.3e   (Hurwitz if < 0)\n', max(real(eJ0)));
fprintf('  min |eig(J0)|  = %.3e\n', min(abs(eJ0)));

M_at_zero = M_func(0);   % = -C_scaled * (J0 \ (B_scaled * diag(W)))

fprintf('norm(M) at w=0:           %.4e\n', norm(M_at_zero));
fprintf('Property 1 (mu) at w=0:   %.4e\n', mu_upper_bound_property1(M_at_zero));
fprintf('Property 2 (gamma) at w=0:%.4e\n', mu_upper_bound_property2(M_at_zero));

T = abs(M_at_zero).^2;
T_pert = T + 0.0005 * ones(q,q);
fprintf('Perron at w=0:            %.4e\n', sqrt(max(abs(eig(T_pert)))));
fprintf('Recursive at w=0:         %.4e\n', mu_upper_bound_recursive(M_at_zero));


fprintf('----------------------------------------\n\n');
fprintf('Computing bounds across %d frequencies...\n', n_grid);
for k = 1:n_grid
    w = w_grid(k);
    s = 1j * w;

    % Build M(w) using the function handle 
    M = M_func(s);

    mu_bound(k)        = mu_upper_bound_property1(M);
    gamma_bound(k)     = mu_upper_bound_property2(M);
    recursive_bound(k) = mu_upper_bound_recursive(M);

    % Perron eigenvalue bound
    T = abs(M).^2;
    epsilon = 0.0005;
    T_pert = T + epsilon * ones(q, q);
    gamma_perron(k) = sqrt(max(abs(eig(T_pert))));

    % MATLAB toolbox mussv 
    [mu_b, ~] = mussv(M, blk);
    mu_mussv(k) = mu_b(1);
end

% Find the maximum  
mu_max        = max(mu_bound);
gamma_max     = max(gamma_bound);
recursive_max = max(recursive_bound);
perron_max    = max(gamma_perron);

fprintf('\n Results \n');
fprintf('Property 1 (mu)     sup = %.4f   ||Delta||_inf < %.4f\n', ...
    mu_max, 1 / mu_max);
fprintf('Property 2 (gamma)  sup = %.4f   ||Delta||_2   < %.4f\n', ...
    gamma_max, 1 / gamma_max);
fprintf('Recursive (gamma)   sup = %.4f   ||Delta||_2   < %.4f\n', ...
    recursive_max, 1 / recursive_max);
fprintf('Perron    (gamma)   sup = %.4f   ||Delta||_2   < %.4f\n', ...
    perron_max, 1 / perron_max);

% Interpret radii in the [-1,1] box 
fprintf('\nBox-uncertainty interpretation (Delta_i in [-1,1]):\n');

inf_radius      = 1 / mu_max;
two_norm_radius = 1 / gamma_max;
sqrt_q          = sqrt(q); 

if inf_radius >= 1
    fprintf('  inf-norm radius = %.4f >= 1.0           : box is covered\n', inf_radius);
else
    fprintf('  inf-norm radius = %.4f < 1.0            : box is NOT covered\n', inf_radius);
end

if two_norm_radius >= sqrt_q
   fprintf('  2-norm radius   = %.4f >= sqrt(q) = %.3f : box is covered\n', two_norm_radius, sqrt_q);
else
    fprintf('  2-norm radius   = %.4f < sqrt(q) = %.3f  : box is NOT covered\n', two_norm_radius, sqrt_q);
end

% Plots 
figure('Color', 'w');
semilogx(w_grid, mu_bound, 'r--', 'LineWidth', 1.5); hold on;
semilogx(w_grid, gamma_bound, 'b-', 'LineWidth', 1.5);
xlabel('\omega'); ylabel('Upper bounds');
title(sprintf('%s: mu vs gamma upper bounds', strrep(model_selection, '_', ' ')));
legend('mu (Property 1)', 'gamma (Property 2)');
grid on; box on;

figure('Color', 'w');
semilogx(w_grid, gamma_bound, 'r--', 'LineWidth', 1.5); hold on;
semilogx(w_grid, gamma_perron, 'b-', 'LineWidth', 1.5);
semilogx(w_grid, recursive_bound, 'k-', 'LineWidth', 1.5);
xlabel('\omega'); ylabel('Upper bounds');
title(sprintf('%s: Gamma upper bounds', strrep(model_selection, '_', ' ')));
legend('(Property 2)', '(Perron bound)','(Recursive bound)');
grid on; box on;


