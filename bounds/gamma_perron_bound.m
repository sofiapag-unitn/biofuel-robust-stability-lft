function gp = gamma_perron_bound(M, blk, epsilon)

% M is partitioned into an nb x nb grid of sub-blocks according to blk 
%       T(i,j) = sigma_bar( M_ij )^2          (= |M_ij|^2 when the block is 1x1)


    n = size(M, 1);

    if nargin < 3 || isempty(epsilon)
        epsilon = 0;
    end
    if nargin < 2 || isempty(blk)
        blk = [-ones(n, 1), zeros(n, 1)];   % default: n independent 1x1 blocks
    end

    % channel index range for each block
    nb  = size(blk, 1);
    rng = cell(nb, 1);
    pos = 1;
    for i = 1:nb
        r      = abs(blk(i, 1));
        rng{i} = pos:(pos + r - 1);
        pos    = pos + r;
    end

    % T(i,j) = sigma_bar( M_ij )^2  over the block partition of M
    T = zeros(nb, nb);
    for i = 1:nb
        for j = 1:nb
            Mij    = M(rng{i}, rng{j});
            T(i,j) = norm(Mij, 2)^2;        % largest singular value squared
        end
    end

    % Perron root of the (regularised) nonnegative matrix
    T  = T + epsilon * ones(nb, nb);
    gp = sqrt(max(abs(eig(T))));
end
