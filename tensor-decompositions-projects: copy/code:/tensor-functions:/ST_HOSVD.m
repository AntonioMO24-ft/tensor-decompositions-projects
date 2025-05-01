function [G, U, V, W, ERR] = ST_HOSVD(X, rank_error)
    % Function to compute the ST-HOSVD or the sequentially truncated HOSVD.
    % This function is usually favored since it fixes two factor matrices 
    % first and then computes the SVD of each one of the other unfoldings 
    % one by one.
    % -----------------------------
    % Inputs: X - tensor we are working with. rank_error - either an array
    % of the desired ranks for the core tensor or a scalar for the desired
    % relative error of the core tensor. 
    % -----------------------------
    % Outputs: G - core tensor of the Tucker decomposition. U,V,W factor
    % matrices of the first, second, and third mode of the tucker
    % decomposition, respectively. ERR - relative error of the decomposition.
    
    [m, n, p] = size(X); % dimensions of X
    chi = norm(X, "fro");  % Frobenius norm of X

    if isscalar(rank_error)
        eps = (rank_error / sqrt(3)) * chi;
        use_error = true;
    else
        q = rank_error(1);
        r = rank_error(2);
        s = rank_error(3);
        use_error = false;
    end
    
    % here we compute the mode-1 unfolding girst and the LLSV of the mode.
    % we fix the other two factor matrices and and U is the LLSV of X_(1)
    X1 = reshape(X, m, []); % Mode-1 unfolding
    if use_error
        [U, ~] = LLSV(X1, eps);
    else
        [U, eps1] = LLSV(X1, q);
    end
    
    % and compress
    Y = ttm(X, U', 1);  % X x_1 U^T

    % same as above now in mode 2
    Y2 = reshape(Y, n, []);
    if use_error
        [V, ~] = LLSV(Y2, eps);
    else
        [V, eps2] = LLSV(Y2, r);
    end

    Z = ttm(Y, V', 2);  % Y x_2 V^T

    Z3 = reshape(Z, p, []); % Mode-3 unfolding
    if use_error
        [W, ~] = LLSV(Z3, eps);
    else
        [W, eps3] = LLSV(Z3, s);
    end

    % now getting G
    G = ttm(Z, W', 3);  % Z ×_3 W^T

    % and computing the error as always
    ERR = sqrt(eps1^2+eps2^2+eps3^2);
end


