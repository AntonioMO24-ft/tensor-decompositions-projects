function [A, B, C, lambda] = CP_ALS3(X, R, max_iter, tol)
    % Function that computes the CP decompositions using the ALS method
    % (alternating least squares). This function only applies for 3-way
    % tensors. See CP_ALSd for general d-way tensors.
    % ------------------
    % Inputs: X - 3-way tensor of mxnxp size, R - desired rank, max_iter -
    % maximum number of iterations, tol - desired tolerance. 
    % ------------------
    % Outputs: A,B,C - factor matrices of the decomposition, lambda - a
    % vector of the weights (can be the norms of the factor matrices).

    % Parameters
    e_old = inf;
    chi = norm(X);
    dimensions = size(X);
    n = dimensions(2); 
    p = dimensions(3);

    % we do not need A or the dimensions of it becuase we work with it
    % first
    B = randn(n, R); % random initializations
    C = randn(p, R);
    S2 = B' * B;
    S3 = C' * C;
    
    % main loop
    for iter = 1:max_iter
        % Updating A
        X1 = double(tenmat(X, 1));  % mode-1 of X
        U1 = X1 * khatrirao(C, B);
        V1 = S2 .* S3;
        A = U1 / V1;
        [A, ~] = ColumnNormalize(A);
        S1 = A' * A;

        % updating B
        X2 = double(tenmat(X, 2));  % mode-2 of X
        U2 = X2 * khatrirao(C, A);
        V2 = S3 .* S1;
        B = U2 / V2;
        [B, ~] = ColumnNormalize(B);
        S2 = B' * B;

        % updating C
        X3 = double(tenmat(X, 3));  % mode-3 of X
        Z = khatrirao(B, A);  % Used in error computation
        U3 = X3 * Z;
        V3 = S2 .* S1;
        C = U3 / V3;
        [C, lambda] = ColumnNormalize(C);
        S3 = C' * C;

        % Computing the error and checking for convergence
        % alpha and beta are flipped from the code
        % lambdas are transposed because of ColumnNormalize
        alpha = lambda * ((S3 .* V3) * lambda'); 
        beta = sum((C' * (X3 * Z)) .* lambda', 'all'); 
        et = sqrt(chi^2 - 2 * alpha + beta);
        et = double(et);

        if iter > 1 && (abs(e_old-et) < tol * chi)
            break;
        end
        e_old = et;
    end
end

% helper function to normalize the columns of the factor matrices
% lambda is a vector that stores the norms of each factor matrix
function [A, lambda] = ColumnNormalize(A)
    [~, r] = size(A);
    lambda = zeros(1, r);

    for j = 1:r
        lambda(j) = norm(A(:, j), 2); 
        if lambda(j) > 0
            A(:, j) = A(:, j) / lambda(j); 
        end
    end
end

