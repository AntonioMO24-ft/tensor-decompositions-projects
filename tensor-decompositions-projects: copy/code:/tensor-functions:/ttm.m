function [Y] = ttm(X,U,k)
    % Function to compute Tensor times matrix in mode k
    % -----------------------------
    % Inputs: X - input tensor. U - input matrix. k - is the mode we want
    % to get the product with respect of X and U, k is in [d].
    % -----------------------------
    % Outputs: Y - resulting tensor.

    X_size = size(X); % Get size of the tensor X
    n_k = X_size(k);  % Dimension of X along mode k
    r = size(U, 1);   % Target dimension after multiplication

    if k == 1 % doing k=1 first is faster and more efficient than including it in the loop below
        % Direct multiplication when k = 1 (mode-1 unfolding)
        X_mat = reshape(X, n_k, []);  % Unfold X along mode-1. kept getting the error so did []
        Y_mat = U * X_mat;  % Multiply with U (r x n_k) * (n_k x M)
        Y = reshape(Y_mat, [r, X_size(2:end)]); % reshape back to tensor form
    else
        % general case for k > 1
        M_k = prod(X_size(1:k-1)); % product of the dimensions before k
        P_k = prod(X_size(k+1:end)); % product of dimensions after k
        X_bar = reshape(X, [M_k, n_k, P_k]); % unfold the tensor into the 3D array

        % here we do matrix multiplication for each of the frontal slices of
        % X
        Y_bar = zeros(M_k, r, P_k); % preallocation
        for l = 1:P_k
            Y_bar(:, :, l) = X_bar(:, :, l) * U'; % multiplying each slice
        end

        % reshaping back to tensor form
        new_size = [X_size(1:k-1), r, X_size(k+1:end)];
        Y = reshape(Y_bar, new_size);
    end
end






