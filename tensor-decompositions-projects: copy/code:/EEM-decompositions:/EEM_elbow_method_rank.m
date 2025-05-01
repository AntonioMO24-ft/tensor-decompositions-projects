% Script that uses the cp-als function from the tensor toolbox to compute 
% the CP model of the EEM tensor for ranks r ∈ { 1, 2, 3, 4, 5 }, using 
% five different initial guesses for each r. Moreover, we also plot the 
% lowest relative error for each rank.

% Parameters
ranks = 1:5; % we will be computing CPD on ranks 1 to 5
initializations = 5; % number of initializations per run
errors = zeros(length(ranks), initializations); % preallocating the errors
best_errors = zeros(length(ranks),1); % storing the best error for all intializations per run
rank_solutions = cell(max(ranks), initializations); % storing the solutions to plot

% Starting the decomposition over each rank
for i = 1:length(ranks)
    r = ranks(i); % current rank we are doing decomposition on
    min_error = inf;
    for j = 1:initializations
        rng(j); % seed for reproducibility
        M = cp_als(X, r);
        % computes an estimate of the best rank-R CP model of a tensor X 
        % using an alternating least-squares algorithm.
        
        % relative error
        X_approx = full(M); % full() reconstructs the tensor
        error = norm(X(:) - X_approx(:)) / norm(X(:)); % computing the error using fro norm
        
        % storing the minimum error over all intializations
        min_error = min(min_error, error);
        errors(i, j) = error; % storing it
        rank_solutions{r, j} = M; 
        % rank solutions is an array containing all runs for every rank to
        % use them in later exercises like 9.12
    end
    
    best_errors(i) = min_error; % the minimum errors
end
save('rank_solutions'); % saving them for use later

% Plot
figure;
plot(ranks, best_errors, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Rank'); ylabel('Lowest Relative Error');
title('CP Decomposition: Lowest Error vs. Rank');
grid on;

% Best rank appears to be r = 3 with lowest error of 0.0331


