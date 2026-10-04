clear;
clc;
close all;
%  Empirical robustness-radius sweep: BOX (inf-norm) and SPHERE (2-norm)
%  Experiments (pick one):
%    'A_total'   Model A, LFT total derivative   m = 6   (kinetic channels)
%    'A_eqparam' Model A, LFT eq.-as-parameter   m = 9   (6 kin + 3 eq coords)
%    'A_bdc'     Model A, BDC ranges             m = 7   (Jacobian entries D)
%    'B_total'   Model B, LFT total derivative   m = 9
%    'B_eqparam' Model B, LFT eq.-as-parameter   m = 12  (9 kin + 3 eq coords; R_bar fixed, not swept)
%    'B_bdc'     Model B, BDC ranges             m = 8

experiment = 'A_bdc';     % <-- SELECT HERE
rho        = 0.20;
alpha_grid = [0.5 0.8 1.0 1.2 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0];
r_grid     = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0];
Nmc        = 500;           % uniform interior draws per alpha
Nsphere    = 500;           % random unit directions per radius
rng(1); % fixed seed
% parse experiment -> model, method, m, certified values
[model, method, m, muInv, gamInv] = parse_experiment(experiment);
sqrt_m = sqrt(m);
% model constants
switch model
case 'A'
        theta_nom = [0.66 0.91 0.14 0.66 0.10 0.50];          % 6 kinetics
        m_kin = 6;  m_eq = 3;  q = 7;
        % Equilibrium box from modelA_openloop_check.m, 6-parameter box
        eqmin = [0.0866 0.2509 0.0469];                       % n,p,b
        eqmax = [0.2694 0.3763 0.2285];
        R_fixed = NaN;   % unused for Model A
case 'B'
        theta_nom = [0.66 0.91 0.14 0.66 0.10 0.50 0.20 100 1.8];  % 9 kinetics
        m_kin = 9;  m_eq = 3;  q = 8;
        % Equilibrium box from modelB_closedloop_check.m, 9-parameter box 
        % R_bar is held FIXED 
        eqmin = [0.0724 0.0957 0.0751];                       % n,p,b
        eqmax = [0.2909 0.2978 0.4954];
        a_R=0.01; k_R=10; g_I=60; I=1000; b_R=2.1;
        R_fixed = (a_R + k_R*(I/(I+g_I)))/b_R;                % closed-form nominal R_bar (~4.4971)
end
mid_eq = (eqmin+eqmax)/2;  hw_eq = (eqmax-eqmin)/2;
% BDC structural matrices and D-box
[Bs, Cs] = bdc_BC(model);
[midD, hwD, Dnom] = bdc_Dbox(model, theta_nom, rho, m_kin);

alpha_phys_kin = 1/rho;
% eq-as-parameter and BDC each hit zero exactly at mid_i/hw_i per channel,
alpha_phys_eq  = min(mid_eq ./ hw_eq);
alpha_phys_bdc = min(midD   ./ hwD);
switch method
case 'kinetic'
    alpha_phys_closed = alpha_phys_kin;
    phys_note = 'raw-parameter wall (1/rho) ONLY -- the equilibrium re-solve is NOT closed-form and typically binds earlier; see the empirical wall below';
    phys_exact = false;
case 'eqparam'
    alpha_phys_closed = min(alpha_phys_kin, alpha_phys_eq);
    phys_note = 'exact: min(kinetic-parameter wall, equilibrium-coordinate wall)';
    phys_exact = true;
case 'bdc'
    alpha_phys_closed = alpha_phys_bdc;
    phys_note = 'exact: min_i(D_i/hw_i), the only physicality gate this method has';
    phys_exact = true;
end
% -------------------------------------------------------------------
% pack constants
S.theta_nom=theta_nom; S.rho=rho; S.model=model;
S.m_kin=m_kin; S.mid_eq=mid_eq; S.hw_eq=hw_eq;
S.midD=midD; S.hwD=hwD; S.B=Bs; S.C=Cs; S.R_fixed=R_fixed;
% nominal anchor
fprintf('==== Experiment %s :  Model %s, method ''%s'',  m = %d ====\n', ...
        experiment, model, method, m);
[mr0, ok0] = evalexp(zeros(1,m), method, S);
fprintf('Anchor (delta=0): max Re eig = %+.4e  (physical=%d)\n', mr0, ok0);
if strcmp(method,'bdc')
    fprintf('Computed D-box (D_min .. D_max):\n');
for i=1:q, fprintf('   D%-2d  [% .5g , % .5g]\n', i, midD(i)-hwD(i), midD(i)+hwD(i)); end
end
fprintf('\n');
% corner set for the box worst-case
if m <= 13
    Sgn = 2*(dec2bin(0:2^m-1) - '0') - 1;   % all 2^m corners
else
    Sgn = sign(randn(8192, m)); Sgn(Sgn==0)=1;
end
nC = size(Sgn,1);
% BOX (inf-norm) sweep
fprintf('--- BOX (inf-norm) sweep : worst over %d corners + %d interior draws ---\n', nC, Nmc);
fprintf('   %7s %14s %10s %10s %12s %10s %14s\n','alpha','worst maxRe','#unst(cor)','#nphys(cor)','P(stab,MC)','#loc.unst','worst(loc)');
Krefine   = 20;    % nearby points tested around each unstable MC hit
eps_local = 0.02;  % size of that local neighborhood, in normalized delta units
alpha_wall_empirical = NaN;   
for a = alpha_grid
    worst=-inf; nun=0; nbad=0;
for c=1:nC
        [mr,ok] = evalexp(a*Sgn(c,:), method, S);
if ~ok, nbad=nbad+1; continue; end
if mr>=0, nun=nun+1; end
if mr>worst, worst=mr; end
end
if isnan(alpha_wall_empirical) && nbad>0     
        alpha_wall_empirical = a;
end
    nok=0; nunm=0;
    n_loc_unst=0; worst_loc=-inf;
for i=1:Nmc
        delta0 = a*(2*rand(1,m)-1);
        [mr,ok] = evalexp(delta0, method, S);
if ~ok, continue; end
        nok=nok+1;
if mr>=0
            nunm=nunm+1;
% found an unstable interior point, probe its neighborhood
for j=1:Krefine
                delta_local = delta0 + eps_local*(2*rand(1,m)-1);
                [mr_loc, ok_loc] = evalexp(delta_local, method, S);
if ok_loc
if mr_loc>=0, n_loc_unst = n_loc_unst+1; end
if mr_loc>worst_loc, worst_loc=mr_loc; end
end
end
end
end
    Pbox=(nok-nunm)/max(nok,1);
    fprintf('   %7.3f %14.4e %10d %10d %12.5f %10d %14.4e\n', a, worst, nun, nbad, Pbox, n_loc_unst, worst_loc);
end
% SPHERE (2-norm) sweep
fprintf('\n--- SPHERE (2-norm) sweep : %d random unit directions per radius ---\n', Nsphere);
fprintf('   %7s %12s %14s %10s %10s\n','r','P(stable)','worst maxRe','#unstable','#nonphys');
for r = r_grid
    worst=-inf; nun=0; nbad=0; nok=0;
for i=1:Nsphere
        u=randn(1,m); u=u/norm(u);
        [mr,ok] = evalexp(r*u, method, S);
if ~ok, nbad=nbad+1; continue; end
        nok=nok+1; if mr>=0, nun=nun+1; end
if mr>worst, worst=mr; end
end
    Psph=(nok-nun)/max(nok,1);
    fprintf('   %7.3f %12.5f %14.3e %10d %10d\n', r, Psph, worst, nun, nbad);
end
% comparison
fprintf('\n--- Certified vs empirical (Model %s, %s, m=%d) ---\n', model, method, m);
fprintf('  BOX    : certified 1/mu_sup    = %.4f   (covered iff >= 1)\n', muInv);
fprintf('  SPHERE : certified 1/gamma_sup = %.4f   (box covered iff >= sqrt(m) = %.4f)\n', gamInv, sqrt_m);
fprintf('  Empirical radius = scale at which P(stable) drops below 1 / nonphysical sets in.\n');

fprintf('\n--- Physicality-aware radius (box, alpha units) ---\n');
if phys_exact
    fprintf('  Closed-form physicality wall = %.4f   (%s)\n', alpha_phys_closed, phys_note);
else
    fprintf('  Closed-form physicality wall = %.4f   (%s)\n', alpha_phys_closed, phys_note);
end
if ~isnan(alpha_wall_empirical)
    fprintf('  Empirical physicality wall (from box sweep #nphys column, grid resolution) = %.4f\n', alpha_wall_empirical);
else
    fprintf('  Empirical physicality wall not reached within tested alpha_grid [%.2f, %.2f]\n', min(alpha_grid), max(alpha_grid));
end
alpha_phys_report = alpha_phys_closed;
if ~phys_exact && ~isnan(alpha_wall_empirical)

    alpha_phys_report = min(alpha_phys_closed, alpha_wall_empirical);
end
radius = min(muInv, alpha_phys_report);
if muInv <= alpha_phys_report
    fprintf('  min(1/mu, alpha_phys) = %.4f  ->  certificate is the binding constraint\n', radius);
else
    fprintf('  min(1/mu, alpha_phys) = %.4f  ->  PHYSICALITY is the binding constraint (1/mu alone overstates the safe radius)\n', radius);
end
% -------------------------------------------------
%% ===================== local functions =====================
function [model, method, m, muInv, gamInv] = parse_experiment(e)
switch e
case 'A_total',   model='A'; method='kinetic'; m=6;  muInv=1.7323; gamInv=3.5962;
case 'A_eqparam', model='A'; method='eqparam'; m=9;  muInv=1.3599; gamInv=1.8641;
case 'A_bdc',     model='A'; method='bdc';     m=7;  muInv=1.5171; gamInv=2.1180;
case 'B_total',   model='B'; method='kinetic'; m=9;  muInv=1.1920; gamInv=3.4245;
case 'B_eqparam', model='B'; method='eqparam'; m=12; muInv=1.0917; gamInv=1.6083;
case 'B_bdc',     model='B'; method='bdc';     m=8;  muInv=0.9616; gamInv=1.8154;
otherwise, error('unknown experiment %s', e);
end
end
function [mr, ok] = evalexp(delta, method, S)
% map a (scaled) normalized perturbation delta in R^m to max Re eig(J).
    mr = NaN; ok = false;
switch method
case 'kinetic'
            theta = S.theta_nom .* (1 + S.rho*delta(:)');
if any(theta<=0), return; end
            [xbar, okx] = solve_eq(theta, S.model);
if ~okx, return; end
            J = Jac(theta, xbar, S.model);
case 'eqparam'
            dk = delta(1:S.m_kin);  de = delta(S.m_kin+1:end);
            theta = S.theta_nom .* (1 + S.rho*dk(:)');
if any(theta<=0), return; end
            xbar_free = S.mid_eq + S.hw_eq.*de(:)';      % equilibrium coords as free channels (n,p,b)
if any(xbar_free<=0), return; end
if strcmp(S.model,'A')
                xbar = xbar_free;                                 % [n p b]
else
                xbar = [xbar_free(1), S.R_fixed, xbar_free(2), xbar_free(3)];  % [n R p b]; R fixed, not swept
end
            J = Jac(theta, xbar, S.model);
case 'bdc'
            D = S.midD + S.hwD.*delta(:)';
if any(D<=0), return; end
            J = S.B*diag(D)*S.C;
end
    mr = max(real(eig(J)));  ok = true;
end
function [xbar, ok] = solve_eq(theta, model)
% re-solve the equilibrium at kinetic vector theta. Returns row vector.
    xbar = []; ok = false;
if strcmp(model,'A')
        a_n=theta(1); d_n=theta(2); g_p=theta(3); b_p=theta(4); a_b=theta(5); d_b=theta(6);
        a_p=0.01; k_p=0.20; g_I=60; I=1000; nmax=1.0;
        p = (a_p + k_p*(I/(I+g_I)))/b_p;
        n = (a_n*g_p/(p+g_p)) / (a_n/nmax + d_n*a_b/(d_b*p));
        b = a_b*n/(d_b*p);
if ~(n>0 && p>0 && b>0 && isfinite(n+p+b)), return; end
        xbar = [n p b];  ok = true;
else
        a_n=theta(1); d_n=theta(2); g_p=theta(3); b_p=theta(4); a_b=theta(5); d_b=theta(6);
        k_p=theta(7); k_b=theta(8); g_R=theta(9);
        a_p=0.01; g_I=60; I=1000; nmax=1.0; a_R=0.01; k_R=10; b_R=2.1;
        R = (a_R + k_R*(I/(I+g_I)))/b_R;
        f = @(x)[ x(1)*(a_n*(1 - x(1)/nmax) - d_n*x(3) - a_n*x(2)/(x(2)+g_p));
                  a_p + k_p/(R/(1+k_b*x(3)) + g_R) - b_p*x(2);
                  a_b*x(1) - d_b*x(2)*x(3) ];

        opts = optimoptions('fsolve','Display','off', ...
'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'StepTolerance',1e-12);
        [xb,fval,flag] = fsolve(f, [0.17;0.17;0.21], opts);
if flag<=0, return; end
if norm(fval) > 1e-8, return; end   % now a real "solver actually failed" gate, not a false-positive trap
        n=xb(1); p=xb(2); b=xb(3);
if ~(n>1e-3 && p>1e-3 && b>1e-3 && isfinite(n+p+b)), return; end
        xbar = [n R p b];  ok = true;
end
end
function J = Jac(theta, xbar, model)
% exact Jacobian (3x3) given kinetics theta and equilibrium coords xbar.
if strcmp(model,'A')
        a_n=theta(1); d_n=theta(2); g_p=theta(3); b_p=theta(4); a_b=theta(5); d_b=theta(6);
        nmax=1.0; n=xbar(1); p=xbar(2); b=xbar(3);
        phi = g_p/(p+g_p)^2;
        J = [ -a_n*n/nmax, -a_n*n*phi, -d_n*n;
               0,          -b_p,        0;
               a_b,        -d_b*b,     -d_b*p ];
else
        a_n=theta(1); d_n=theta(2); g_p=theta(3); b_p=theta(4); a_b=theta(5); d_b=theta(6);
        k_p=theta(7); k_b=theta(8); g_R=theta(9);
        nmax=1.0; n=xbar(1); R=xbar(2); p=xbar(3); b=xbar(4);
        phi = g_p/(p+g_p)^2;
        denomR = R/(1+k_b*b) + g_R;
        psi_b  = k_p*R*k_b/(denomR^2*(1+k_b*b)^2);
        J = [ -a_n*n/nmax, -a_n*n*phi, -d_n*n;
               0,          -b_p,        psi_b;
               a_b,        -d_b*b,     -d_b*p ];
end
end
function D = Dvec(theta, xbar, model)
% BDC diagonal channels from kinetics and equilibrium.
if strcmp(model,'A')
        a_n=theta(1); d_n=theta(2); g_p=theta(3); b_p=theta(4); a_b=theta(5); d_b=theta(6);
        nmax=1.0; n=xbar(1); p=xbar(2); b=xbar(3);
        phi = g_p/(p+g_p)^2;
        D = [ a_n*n/nmax, a_n*n*phi, d_n*n, b_p, a_b, d_b*b, d_b*p ];
else
        a_n=theta(1); d_n=theta(2); g_p=theta(3); b_p=theta(4); a_b=theta(5); d_b=theta(6);
        k_p=theta(7); k_b=theta(8); g_R=theta(9);
        nmax=1.0; n=xbar(1); R=xbar(2); p=xbar(3); b=xbar(4);
        phi = g_p/(p+g_p)^2;
        denomR = R/(1+k_b*b) + g_R;
        psi_b  = k_p*R*k_b/(denomR^2*(1+k_b*b)^2);
        D = [ a_n*n/nmax, a_n*n*phi, d_n*n, b_p, psi_b, a_b, d_b*b, d_b*p ];
end
end
function [Bs, Cs] = bdc_BC(model)
% structural matrices so that J = Bs*diag(D)*Cs.
if strcmp(model,'A')
        Bs = zeros(3,7);
        Bs(1,1)=-1; Bs(1,2)=-1; Bs(1,3)=-1; Bs(2,4)=-1; Bs(3,5)=1; Bs(3,6)=-1; Bs(3,7)=-1;
        Cs = zeros(7,3);
        Cs(1,1)=1; Cs(2,2)=1; Cs(3,3)=1; Cs(4,2)=1; Cs(5,1)=1; Cs(6,2)=1; Cs(7,3)=1;
else
        Bs = zeros(3,8);
        Bs(1,1)=-1; Bs(1,2)=-1; Bs(1,3)=-1; Bs(2,4)=-1; Bs(2,5)=1; Bs(3,6)=1; Bs(3,7)=-1; Bs(3,8)=-1;
        Cs = zeros(8,3);
        Cs(1,1)=1; Cs(2,2)=1; Cs(3,3)=1; Cs(4,2)=1; Cs(5,3)=1; Cs(6,1)=1; Cs(7,2)=1; Cs(8,3)=1;
end
end
function [midD, hwD, Dnom] = bdc_Dbox(model, theta_nom, rho, m_kin)
% Step 1: list every combination of +1/-1 for m_kin parameters
% There are 2^m_kin combinations
    numCombos = 2^m_kin;
    signCombos = zeros(numCombos, m_kin);
for c = 1:numCombos
        code = c - 1;              % combination index, counting from 0
for k = 1:m_kin
            bitValue = mod(code, 2);   % 0 or 1
            code = floor(code / 2);    % move to the next parameter's bit
if bitValue == 0
                signCombos(c, k) = -1;  % this parameter goes down
else
                signCombos(c, k) = 1;   % this parameter goes up
end
end
end
% Step 2: compute D at the nominal (unperturbed) parameters
    [xb0, ~] = solve_eq(theta_nom, model);
    Dnom = Dvec(theta_nom, xb0, model);
% Step 3: set up trackers for the smallest/largest D seen
    Dmin = inf(size(Dnom));
    Dmax = -inf(size(Dnom));
% Step 4: try every combination
for c = 1:numCombos
% perturb each kinetic parameter up or down by rho (e.g. 20%)
        theta = theta_nom .* (1 + rho * signCombos(c, :));
% skip this combination if it made a parameter zero or negative
if any(theta <= 0)
continue;
end
% re-solve the equilibrium for this perturbed parameter set
        [xb, ok] = solve_eq(theta, model);
if ~ok
continue;   % skip if the solver failed to converge
end
% compute D at this perturbed point
        D = Dvec(theta, xb, model);
% update the running min/max, one entry at a time
for k = 1:length(D)
if D(k) < Dmin(k)
                Dmin(k) = D(k);
end
if D(k) > Dmax(k)
                Dmax(k) = D(k);
end
end
end
% Step 5: convert min/max range into midpoint + half-width
    midD = (Dmin + Dmax) / 2;
    hwD  = (Dmax - Dmin) / 2;
end
