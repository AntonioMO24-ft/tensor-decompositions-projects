function [G, U, V, W, ERR] = HOSVD3(X, rank_error)
    % HOSVD3 computes the Tucker decomposition of a 3-way tensor X. Same as
    % HOSVD function but this one is specifically for a 3-way tensor.
    % -----------------------------
    % Inputs: X - tensor we are working with. rank_error - either an array
    % of the desired ranks for the core tensor or a scalar for the desired
    % relative error of the core tensor. 
    % -----------------------------
    % Outputs: G - core tensor of the Tucker decomposition. U,V,W factor
    % matrices of the first, second, and third mode of the tucker
    % decomposition, respectively. ERR - relative error of the decomposition. 
    
    [m, n, p] = size(X); % getting the dimensions of the tensor X
    chi = norm(X, "fro"); % computing the Frobenius norm 
    
    % checking if an error tolerance was given
    if isscalar(rank_error)
        eps = (rank_error / sqrt(3)) * chi; % adjusting the error for the threshold
        use_error = true;
    else % if an array of the dimensions of the core tensors are given
        q = rank_error(1);
        r = rank_error(2);
        s = rank_error(3);
        use_error = false;
    end
    
    % here we compute the unfoldings of the tensor
    X1 = reshape(X, m, n * p); % mode-1 unfolding
    X2 = reshape(X, n, m * p); % mode-2 unfolding
    X3 = reshape(X, p, m * n); % mode-3 unfolding

    % here we compute the LLSV's and allocated them for each factor matrix
    if use_error
        % Use error tolerance to compute rank
        [U, ~] = LLSV(X1, eps);
        [V, ~] = LLSV(X2, eps);
        [W, ~] = LLSV(X3, eps);
    else
        % Use specific rank for each mode
        [U, ~] = LLSV(X1, q);  % Rank for mode-1 
        [V, ~] = LLSV(X2, r);  % Rank for mode-2
        [W, ~] = LLSV(X3, s);  % Rank for mode-3
    end
    
    % here we perform the ttm's to compute the core tensor
    G1 = ttm(X, U', 1); % Mode-1 multiplication
    G2 = ttm(G1, V', 2); % Mode-2 multiplication
    G = ttm(G2, W', 3); % Mode-3 multiplication

    % and compute the error ERR
    ERR = sqrt(chi^2 - (norm(G, "fro"))^2);
end

