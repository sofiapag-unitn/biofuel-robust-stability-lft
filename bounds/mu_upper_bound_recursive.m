function gamma = mu_upper_bound_recursive(M, blk)
%
%   [-r 0] = r-times repeated REAL    scalar,
%   [ r 0] = r-times repeated COMPLEX scalar.

% For a block of size r the per-block step-2/step-4 LMIs use:
%   H_i : n x r selector for the block's r channels,
%   M_i : r x n  (the block's rows of M),
%   d_i : r x r  Hermitian, d_i >= I_r           
%   g_i : r x r  Hermitian (real block) or 0     
% and step 5 rescales M_k by the block matrix square root sqrtm(d_i).

    q = size(M, 1);
    n = size(M, 1);

    if nargin < 2 || isempty(blk)
        blk = [-ones(q, 1), zeros(q, 1)];   % default: q independent real scalars
    end

    [rng, is_real_blk] = block_ranges(blk);
    nb = numel(rng);

    s = norm(M);
    if s == 0, gamma = 0; return; end
    Mn  = M / s;
    M_k = Mn;

    max_iters      = 8;
    alpha_hat_min  = inf;
    alpha_hat_prev = inf;

    for k = 1:max_iters

        % STEP 2: bisection of alpha_i for each block
        alpha_vec = zeros(nb, 1);
        for i = 1:nb
            alpha_vec(i) = solve_step2_block_bisection(M_k, rng{i}, is_real_blk(i), n);
        end

        % STEP 3
        alpha_hat = max(alpha_vec);
        if alpha_hat < alpha_hat_min
            alpha_hat_min = alpha_hat;
        end
        if alpha_hat_min < 1e-5 || abs(alpha_hat_prev - alpha_hat) < 1e-4
            break;
        end
        alpha_hat_prev = alpha_hat;

        % STEP 4: maximize the block scaling D_i (block-by-block)
        alpha_target = (alpha_vec + alpha_hat) / 2;
        D_sqrt = eye(q);
        for i = 1:nb
            Di = solve_step4_block(M_k, alpha_target(i), rng{i}, is_real_blk(i), n);
            D_sqrt(rng{i}, rng{i}) = sqrtm(Di);    % block matrix square root
        end

        % STEP 5: rescale M
        M_k = D_sqrt * M_k / D_sqrt;

    end

    gamma = s * sqrt(max(0, alpha_hat_min));
end

%--------------------------------------------------------------------------
function [rng, is_real_blk] = block_ranges(blk)
% Channel index range and real/complex flag for each block in blk.
    nb          = size(blk, 1);
    rng         = cell(nb, 1);
    is_real_blk = false(nb, 1);
    pos = 1;
    for i = 1:nb
        r              = abs(blk(i, 1));
        rng{i}         = pos:(pos + r - 1);
        is_real_blk(i) = (blk(i, 1) < 0);   % negative entry -> real block
        pos            = pos + r;
    end
end

%--------------------------------------------------------------------------
function alpha_opt = solve_step2_block_bisection(M_k, idx, is_real, n)
    step2 = build_step2_checker(M_k, idx, is_real, n);
    alpha_lo = 0;
    alpha_hi = 10;
    while ~step2_feasible(step2, alpha_hi)
        alpha_hi = alpha_hi * 10;
        if alpha_hi > 1e6
            alpha_opt = inf;
            return;
        end
    end
    tol = 1e-4;
    while (alpha_hi - alpha_lo) > tol
        alpha_mid = 0.5 * (alpha_lo + alpha_hi);
        if step2_feasible(step2, alpha_mid)
            alpha_hi = alpha_mid;
        else
            alpha_lo = alpha_mid;
        end
    end
    alpha_opt = alpha_hi;
end

%--------------------------------------------------------------------------
function step2 = build_step2_checker(M_k, idx, is_real, n)
    yalmip('clear');
    alpha = sdpvar(1, 1);
    r  = numel(idx);

    Hi = zeros(n, r);                 % n x r selector for this block
    for c = 1:r, Hi(idx(c), c) = 1; end
    Mi = M_k(idx, :);                 % r x n  (block rows of M)
    Li = [Mi', Hi];                   % n x 2r

    d_i = sdpvar(r, r, 'hermitian', 'complex');   % Set_D : Hermitian, d_i >= I_r
    if is_real
        g_i = sdpvar(r, r, 'hermitian', 'complex');   % Set_G : Hermitian (real block)
    else
        g_i = zeros(r, r);                            % complex block -> G = 0
    end

    Term1 = Mi'*d_i*Mi + 1i*(Mi'*g_i*Hi' - Hi*g_i*Mi);
    Term2 = alpha * (Hi*(d_i - eye(r))*Hi' + eye(n));
    Qi = Li' * (Term1 - Term2) * Li;
    Qi = (Qi + Qi') / 2;

    Constraints = [Qi <= 0, d_i >= eye(r)];
    
    opts = sdpsettings('solver','sedumi','verbose',0);
    
    step2 = optimizer(Constraints, [], opts, alpha, d_i);
end

%--------------------------------------------------------------------------
function tf = step2_feasible(step2, alpha_val)
    [~, errorcode] = step2(alpha_val);
    tf = (errorcode == 0);
end

%--------------------------------------------------------------------------
function Di = solve_step4_block(M_k, alpha_tgt, idx, is_real, n)
    yalmip('clear');
    r  = numel(idx);

    Hi = zeros(n, r);
    for c = 1:r, Hi(idx(c), c) = 1; end
    Mi = M_k(idx, :);
    Li = [Mi', Hi];

    d_i = sdpvar(r, r, 'hermitian', 'complex');
    if is_real
        g_i = sdpvar(r, r, 'hermitian', 'complex');
    else
        g_i = zeros(r, r);
    end
    tau = sdpvar(1, 1);

    Term1 = Mi'*d_i*Mi + 1i*(Mi'*g_i*Hi' - Hi*g_i*Mi);
    Term2 = alpha_tgt * (Hi*(d_i - eye(r))*Hi' + eye(n));
    Qi = Li' * (Term1 - Term2) * Li;
    Qi = (Qi + Qi') / 2;

    % Property 4, step 4: maximize tau_i subject to H_i^H D_i H_i >= tau_i I
    Constraints = [Qi <= 0, d_i >= tau*eye(r), tau >= 1, d_i <= 100*eye(r)];
    
    opts = sdpsettings('solver','sedumi','verbose',0);
    
    sol  = optimize(Constraints, -tau, opts);          % maximize tau_i
    if sol.problem == 0 || sol.problem == 3 || sol.problem == 4
        Di = value(d_i);
    else
        Di = eye(r);
    end
end
