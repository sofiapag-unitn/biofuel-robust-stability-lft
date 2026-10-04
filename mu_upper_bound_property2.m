function gamma = mu_upper_bound_property2(M, blk)
% Property 2 (Alamo & Dormido 2001) upper bound to gamma_Delta(M).
%
% Robust rewrite of the original:
%   - NO YALMIP 'optimizer' object (that precompiled path is what crashes
%     SeDuMi on complex/Hermitian models -> "Unrecognized variable y_s").
%   - Each feasibility test is a fresh optimize() call.
%   - Every complex-Hermitian PSD constraint is REAL-EMBEDDED, so SeDuMi
%     only ever sees a real symmetric cone:
%         H >= 0   <=>   [Re(H) -Im(H); Im(H) Re(H)] >= 0
%   - Solver fallback + try/catch so one bad solve can't abort the sweep.
%
%   blk(k,:) = [-r 0]  repeated REAL scalar block of size r
%   blk(k,:) = [ r 0]  repeated COMPLEX scalar block of size r

n = size(M, 1);
if nargin < 2 || isempty(blk)
    blk = [-ones(n,1), zeros(n,1)];
end

s = norm(M);
if s == 0
    gamma = 0;
    return;
end
Mn = M / s;

alpha_lo = 0;
alpha_hi = 10;
tol      = 1e-4;

% expand the upper bracket until feasible
while ~is_feasible_prop2(Mn, blk, alpha_hi)
    alpha_hi = 10 * alpha_hi;
    if alpha_hi > 1e8
        warning('Property 2: LMI infeasible up to alpha = %g.', alpha_hi);
        gamma = NaN;
        return;
    end
end

% bisection
while (alpha_hi - alpha_lo) > tol
    alpha_mid = 0.5 * (alpha_lo + alpha_hi);
    if is_feasible_prop2(Mn, blk, alpha_mid)
        alpha_hi = alpha_mid;
    else
        alpha_lo = alpha_mid;
    end
end

gamma = s * sqrt(max(0, alpha_hi));
end
% =======================================================================
function tf = is_feasible_prop2(M, blk, alpha_val)
% Build and solve the Property-2 LMI feasibility problem for a fixed alpha.
% Returns true iff feasible.

n  = size(M, 1);
nb = size(blk, 1);

% channel range occupied by each block
rng = cell(nb, 1);
pos = 1;
for b = 1:nb
    rb = abs(blk(b, 1));
    rng{b} = pos:(pos + rb - 1);
    pos = pos + rb;
end

% shared D in D_Delta (block-diagonal Hermitian), normalized D >= I
Db = cell(nb, 1);
for b = 1:nb
    rb = abs(blk(b, 1));
    Db{b} = sdpvar(rb, rb, 'hermitian', 'complex');
end
D = blkdiag(Db{:});

Constraints = herm2real(D - eye(n)) >= 0;          % D >= I

for i = 1:nb
    ri = abs(blk(i, 1));

    % Gamma_i = Pii * D  (D restricted to block i)
    Pii = zeros(n, n);
    Pii(rng{i}, rng{i}) = eye(ri);
    Gamma_i = Pii * D;

    % E_i in D_Delta (block-diagonal Hermitian PSD)
    Eb = cell(nb, 1);
    for b = 1:nb
        rb = abs(blk(b, 1));
        Eb{b} = sdpvar(rb, rb, 'hermitian', 'complex');
    end
    Ei = blkdiag(Eb{:});

    % F_i in G_Delta (Hermitian on repeated-real blocks, zero otherwise)
    Fb = cell(nb, 1);
    for b = 1:nb
        rb = abs(blk(b, 1));
        if blk(b, 1) < 0
            Fb{b} = sdpvar(rb, rb, 'hermitian', 'complex');
        else
            Fb{b} = zeros(rb, rb);
        end
    end
    Fi = blkdiag(Fb{:});

    LMI = alpha_val*(D + Ei) - M'*(Gamma_i + Ei)*M - 1i*(M'*Fi - Fi*M);

    Constraints = [ Constraints, ...
                herm2real(LMI) >= 1e-5*eye(2*n), ...
                herm2real(Ei)  >= 0 ];
end

tf = solve_feasibility(Constraints, alpha_val);
end
% =======================================================================
function R = herm2real(H)
% Real embedding of a Hermitian matrix: H >= 0  <=>  R >= 0.
R = [real(H), -imag(H); imag(H), real(H)];
R = 0.5*(R + R.');          % symmetrize to kill round-off asymmetry
end
% =======================================================================
function tf = solve_feasibility(Constraints, alpha_val)
% Try a sequence of SDP solvers; problem == 0 means feasible.
solver_list = {'sedumi', 'sdpt3', 'mosek', 'scs'};
tf = false;
for k = 1:numel(solver_list)
    opts = sdpsettings('verbose', 0, 'solver', solver_list{k}, 'cachesolvers', 1);
    try
        d = optimize(Constraints, [], opts);
    catch ME
        fprintf('    [prop2 debug] alpha=%.6g  solver=%s threw: %s\n', ...
                alpha_val, solver_list{k}, ME.message);
        continue;       % solver missing or crashed -> next one
    end
    if d.problem == -3 || d.problem == -4
        continue;       % solver not found/available -> next one
    end
    if d.problem ~= 0
        fprintf('    [prop2 debug] alpha=%.6g  solver=%s  d.problem=%d (%s)\n', ...
                alpha_val, solver_list{k}, d.problem, yalmiperror(d.problem));
    end
    tf = (d.problem == 0);
    return;
end
warning('Property 2: no SDP solver succeeded on this instance.');
end