% Script to find the best lambda values of the best rank-3 tensor computed
% in EEM_elbow_method_rank.m script. Rank-3 seems to be the best rank.
% Lambdas are the weights of each of the factor matrices.

num_inits = 5; % Number of different initializations
best_error = inf;
best_M = []; % storing the best rank-3 factorization

for j = 1:num_inits
    rng(j); % seed
    M = cp_als(X, 3);
    
    % relative error
    X_approx = full(M);
    error = norm(X(:) - X_approx(:)) / norm(X(:));
    
    % here we store the best error and the best model for a rank 3
    if error < best_error
        best_error = error;
        best_M = M;
    end
end

% .lambda extracts the lambda values for best_M
lambda_values = best_M.lambda;
disp(lambda_values);


% Output:
% Lambda values for rank-3 CP decomposition:
%    1.0e+07 *
% 
%     1.7774
%     0.9509
%     0.4266
