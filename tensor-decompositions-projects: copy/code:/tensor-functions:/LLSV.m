function [W, ERR] = LLSV(Y, rank_error)
    % LLSV computes function that outputs the leading left singular vectors
    % of Y and its relative error, i.e the SVD of the given matrix. 
    % -----------------------------
    % Inputs: Y - input matrix which is the third mode unfolding of tensor 
    % X times V times U. rank_error - either an array of the desired ranks 
    % for the core tensor or a scalar for the desired relative error of the 
    % core tensor. 
    % -----------------------------
    % Outputs: W - solution matrix to the LLSV problem. ERR - relative
    % error.
    
    [U, S, ~] = svd(Y, 'econ');  % here we define the economy SVD of Y
    sigma = diag(S);  % take the singular values from S
    
    % Determine rank r
    if isscalar(rank_error) % if we are given an error tolerance
        eps = rank_error;
        squared_sum = cumsum(sigma(end:-1:1).^2); % cumulative sum computation
        r = find(squared_sum <= eps^2, 1, 'last'); % the smallest r satisfying the condition so the minimum
    else
        r = rank_error; % If we are given a rank
    end

    % now W is the mxr orthogonal matrix that is the solution to the
    % optimization problem
    W = U(:, 1:r); % getting the r leading left singular vectors of matrix Y
    ERR = sqrt(sum(sigma(r+1:end).^2)); % computing the error term
end
