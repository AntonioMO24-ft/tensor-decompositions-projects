% Script that with 50 random initializations, computes the rank-3 CP
% decompositions of the EEM tensor. It also outputs the range of the
% relative errors for the 50 initializations and how often we had a
% relative error of 10% or lower.

% parameters
rank = 3;         
num_trials = 50;  
errors = zeros(num_trials, 1);  % Storing errors
rank_3_solutions = cell(num_trials, 1); % storing decompositions

for trial = 1:num_trials
    rng(trial); % seed
    [M, ~, output] = cp_als(X, rank); % computing cp decomposition
    
    % relative rror
    X_approx = full(M);
    error = norm(X(:) - X_approx(:)) / norm(X(:));
    errors(trial) = error;
    
    % storing the decomposition
    rank_3_solutions{trial} = M;
end

min_error = min(errors); % minimum error
max_error = max(errors); % maximum error
error_range = max_error - min_error; % range
cutoff = sum(errors <= min_error + 0.01); % how many have relative error of 0.01

% results
fprintf('Range of relative errors: [%f, %f]\n', min_error, max_error);
fprintf('Number of trials within 0.01 of the lowest error: %d out of %d\n', cutoff, num_trials);

% Output:
% Range of relative errors: [0.033048, 0.136406]
% Number of trials within 0.01 of the lowest error: 33 out of 50
