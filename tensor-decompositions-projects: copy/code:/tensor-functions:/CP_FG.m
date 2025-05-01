function [f, g] = CP_FG(X, chi, v, r)
    % Function that calculates the function value and gradient value for a
    % given tensor. After obtaining f and g using this function we can run
    % an optimization method to minimize the function value. This is
    % CP-OPT.
    % -------------------
    % Inputs: X - data tensor of size mxnxp, chi - precomputed frobenius
    % norm of the tensor, v - vector of the vectorized forms of the factor
    % matrices stacked, r - desired rank.
    % -------------------
    % Outputs: f - function value, g - gradient vector
    
    dimensions = size(X); % getting the dimensions of the tensor
    m = dimensions(1);
    n = dimensions(2);
    p = dimensions(3);

    % here we get the factor matrices from the given vector v
    A = reshape(v(1:m*r), m, r);
    B = reshape(v(m*r+1 : m*r+n*r), n, r);
    C = reshape(v(m*r+n*r+1:end), p, r);

    % following the algo, we compute the respective gram matrices
    S1 = A' * A;
    S2 = B' * B;
    S3 = C' * C;

    % and the gradients
    G1 = A * (S3 .* S2) - double(tenmat(X, 1)) * khatrirao(C, B);
    G2 = B * (S3 .* S1) - double(tenmat(X, 2)) * khatrirao(C, A);

    % here we compute the hadamard product and the KRP
    V3 = S2 .* S1;
    U3 = double(tenmat(X, 3)) * khatrirao(B, A);

    % computing the final gradient
    G3 = C * V3 - U3;

    % computing the objective function
    f = 0.5 * chi - sum(sum(C .* U3)) + 0.5 * sum(sum(S3 .* V3));

    % vectorizing the gradients to return the vector g
    g = [G1(:); G2(:); G3(:)];
end
