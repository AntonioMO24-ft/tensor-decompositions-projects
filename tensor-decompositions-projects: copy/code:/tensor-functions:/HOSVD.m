function [G, U, ERR] = HOSVD(X, rank_error)
    % Function that computes the HOSVD for d-way tensors. This function
    % should not be chosen over ST-HOSVD since the sequentially truncated
    % version is cheaper and faster. Here rank_error depends on if we want
    % a core tensor with a desired relative error or a desired rank.
    % -----------------------------
    % Inputs: X - tensor we are working with. rank_error - either an array
    % of the desired ranks for the core tensor or a scalar for the desired
    % relative error of the core tensor. 
    % -----------------------------
    % Outputs: G - core tensor of the Tucker decomposition. U - cell with
    % the d factor matrices of the decomposition. ERR - relative error of
    % the decomposition.
    
    dimensions = size(X); % getting the dimensions of the tensor X
    d = length(dimensions); % now we need the dimensions for a d-way tensor
    chi = norm(X, "fro"); % computing the Frobenius norm 
    % checking if an error tolerance was given
    if isscalar(rank_error)
        eps = (rank_error / sqrt(d)) * chi; 
        use_error = true;
    else
        r = rank_error; % if a rank was given it a vector with the dimensions of the rank
        use_error = false;
    end

    U = cell(1, d); % array to store the factor matrices
    errors = zeros(1, d); % array to store the errors of each factor matrix

    for k = 1:d % now we compute the SVD for each mode
        % using reshape and permute to get the k-mode unfolding
        Xk = reshape(permute(X, [k, 1:k-1, k+1:d]), dimensions(k), []);
        if use_error
            [U{k}, errors(k)] = LLSV(Xk, eps); % for each of the factor matrices
        else
            [U{k}, errors(k)] = LLSV(Xk, r(k));
        end
    end
    
    % here we compute the core tensor 
    G = X; % starting G as the tensor to start the loop
    for k = 1:d
        G = ttm(G, U{k}', k); % at each iteration we are doing the k-wise ttm's
    end
    % and computing the approximation error as before
    ERR = sqrt(chi^2 - (norm(G, "fro"))^2);
end
