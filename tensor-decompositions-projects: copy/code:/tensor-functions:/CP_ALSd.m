function [A_factors, lambda] = CP_ALSd(X, R, tol, max_iters)
    % Function that computes the CP decompositions using the ALS method
    % (alternating least squares).
    % ------------------
    % Inputs: X - d-way tensor of mxnxp size, R - desired rank, max_iter -
    % maximum number of iterations, tol - desired tolerance. 
    % ------------------
    % Outputs: A_factors - cell containing the d factor matrices of the 
    % decomposition, lambda - a vector of the weights (can be the norms of 
    % the factor matrices).

    % Parameters
    e_old = inf;
    d = ndims(X);
    dimensions = size(X);
    chi = norm(X);
    
    % Preallocation
    A_factors = cell(d, 1);
    S = cell(d, 1);

    % Initializing the Sk's
    for k = 2:d
        A_factors{k} = randn(dimensions(k), R);
        % A_factors{k} = normalize_columns(A_factors{k});
        S{k} = A_factors{k}' * A_factors{k};
    end
    
    % Main loop
    for iter = 1:max_iters
        for k = 1:d
            % KRP except for the k-mode unfolding
            kr_factors = A_factors([1:k-1, k+1:end]);
            Z = khatrirao(kr_factors(end:-1:1));  % reversed

            Xk = double(tenmat(X, k)); % mode-k unfolding

            % loop to define the Vk's
            Vk = eye(R);
            for j = [1:k-1, k+1:d]
                Vk = Vk .* S{j};
            end
            A_factors{k} = Xk * Z / Vk; % solver

            % normalizing columns
            [A_factors{k}, lambda] = ColumnNormalize(A_factors{k});
            S{k} = A_factors{k}' * A_factors{k};
        end

        % Preallocation to compute the error
        Xd = double(tenmat(X, d));
        Z = khatrirao(A_factors{1:d-1});
        Ud = Xd * Z;
        Vd = eye(R);
        for j = 1:d-1
            Vd = Vd .* S{j};
        end


        alpha = lambda * ((S{d} .* Vd) * lambda');
        beta = sum((A_factors{d}' * Ud) .* lambda', 'all');
        et = sqrt(chi^2 - 2 * alpha + beta);

        if iter > 1 && abs(e_old - et) < tol * chi
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

