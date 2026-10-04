function mu = mu_upper_bound_property1(M, blk)
% Property 1 (Fan et al. upper bound; Alamo & Dormido 2001, Property 1):
%   if there exist alpha in R, D in D_Delta, G in G_Delta with
%       M^H D M + j(M^H G - G M) < alpha * D
%   then  mu_Delta(M) <= sqrt(max(0, alpha)).

%   blk(k,:) = [-r 0]  -> repeated REAL scalar block of size r
%                         (D block: r x r Hermitian PSD ;  G block: r x r Hermitian)
%   blk(k,:) = [ r 0]  -> repeated COMPLEX scalar block of size r
%                         (D block: r x r Hermitian PSD ;  G block: 0)

% If blk is omitted, defaults to n single 1x1 real scalar blocks (the original
% behaviour: D = diag(d), G = diag(g))

n = size(M, 1);
if nargin < 2 || isempty(blk)
    blk = [-ones(n,1), zeros(n,1)];   % all 1x1 real scalar blocks
end

s = norm(M);
if s == 0, mu = 0; return; end
Mn = M / s;

feas = build_checker_p1(Mn, blk);

alpha_lo = 0;
alpha_hi = 10;
tol = 1e-4;
while ~is_feasible(feas, alpha_hi)
    alpha_hi = 10 * alpha_hi;
    if alpha_hi > 1e6
        warning('Property 1 LMI infeasible up to alpha = %g.', alpha_hi);
        mu = NaN;
        return;
    end
end
while (alpha_hi - alpha_lo) > tol
    alpha_mid = 0.5 * (alpha_lo + alpha_hi);
    if is_feasible(feas, alpha_mid)
        alpha_hi = alpha_mid;
    else
        alpha_lo = alpha_mid;
    end
end
mu = s * sqrt(max(0, alpha_hi));
end
% -----------------------------------------------------------------------
function feas = build_checker_p1(M, blk)
    n  = size(M, 1);
    nb = size(blk, 1);
    alpha = sdpvar(1, 1);

    % Build block-diagonal D and G with one Hermitian block per uncertainty
    % block
    Dblocks = cell(nb, 1);
    Gblocks = cell(nb, 1);
    for b = 1:nb
        r = abs(blk(b, 1));                          % block size
        Dblocks{b} = sdpvar(r, r, 'hermitian', 'complex');   
        if blk(b, 1) < 0                             % repeated REAL block -> G Hermitian
            Gblocks{b} = sdpvar(r, r, 'hermitian', 'complex');
        else                                         % complex block -> G = 0
            Gblocks{b} = zeros(r, r);
        end
    end
    D = blkdiag(Dblocks{:}); %assembles the block-diagonal
    G = blkdiag(Gblocks{:});

    LMI = alpha*D - M'*D*M - 1i*(M'*G - G*M);
    Constraints = [ D >= 1e-6*eye(n), LMI >= 1e-6*eye(n) ];
 
    opts = sdpsettings('solver','sedumi','verbose',0);
    
    feas = optimizer(Constraints, [], opts, alpha, D);
end
% -----------------------------------------------------------------------
function tf = is_feasible(feas, alpha_val)
    [~, errorcode] = feas(alpha_val);
    tf = (errorcode == 0);
end