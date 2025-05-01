function [G, U, V, W, ERR, num_iters] = HOOI(X, U, V, W, tau, MAXITERS)
    % Function to compute the Higher-Order Orthogonal Iteration for 3-way
    % Tensors. This algorithm is an iterative algorithm compared to HOSVD
    % and ST-HOSVD which are not. 
    % -----------------------------
    % Inputs: X - tensor we are working with. U,V,W - the left-leading
    % singular vectors of the mode-1 unfolding of X, mode-2 unfolding of X,
    % mode-3 unfolding of X, respectively. tau - tolerance. MAXITERS -
    % maximum amount of iterations.
    % -----------------------------
    % Outputs: G - core tensor of the Tucker decomposition of X. U,V,W -
    % factor matrices of X. ERR - final relative error. num-iters - number
    % of iterations until convergence.
    
    chi = norm(X, "fro");  
    ERR = inf;  
    num_iters = 0; 

    for t = 1:MAXITERS
        num_iters = t;  % Track iterations
        
        % Mode-1 update
        Y1 = ttm(ttm(X, V', 2), W', 3);
        [U, ~] = LLSV(reshape(Y1, size(Y1,1), []), size(U,2));

        % Mode-2 update
        Y2 = ttm(ttm(X, U', 1), W', 3);
        [V, ~] = LLSV(reshape(Y2, size(Y2,2), []), size(V,2));

        % Mode-3 update
        Y3 = ttm(ttm(X, U', 1), V', 2);
        [W, ~] = LLSV(reshape(Y3, size(Y3,3), []), size(W,2));

        % Compute core tensor G
        G = ttm(Y3, W', 3);
        
        % Compute error
        ERR_t = sqrt(chi^2 - norm(G(:))^2);

        % **Ensure at least 5 iterations before stopping**
        if t > 5 && abs(ERR_t - ERR) < tau * chi
            break;
        end
        ERR = ERR_t;
    end
end
