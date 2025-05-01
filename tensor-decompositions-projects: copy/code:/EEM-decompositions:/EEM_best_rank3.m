% Script to compute the best rank-3 CP decomposition of the EEM tensor.
% This uses CP_ALS3.

clear
load('EEM18.mat')

rank = 3; % CP rank
num_trials = 30;

errors = zeros(num_trials, 1);       % Store relative errors
solutions = cell(num_trials, 1);     % Store CP decompositions

for i = 1:num_trials
    rng(i); 
    [A, B, C, lambda] = CP_ALS3(X, rank, 50, 1e-4);

    % Create ktensor from factors and weights
    M = ktensor(lambda(:), {A, B, C});

    
    M = arrange(M);
    % Fix the signs. The key to having nonnegative are these two lines
    
    M = fixsigns(M);

    % Compute approximation error
    X_approx = full(M);  % reconstruct full tensor
    errors(i) = norm(X(:) - X_approx(:)) / norm(X(:));

    % Store CP model
    solutions{i} = M;
end

% Find best CP model
[best_error, best_idx] = min(errors);
best_model = solutions{best_idx};

% Extract factors
A = best_model.U{1}; B = best_model.U{2}; C = best_model.U{3};
R = size(A, 2);

% Plot factor matrices
figure('Name', 'Best CP Decomposition');
for r = 1:R
    subplot(3, R, r);
    bar(A(:, r));
    title(sprintf('Sample %d', r));
    ylabel('Weight');

    subplot(3, R, R + r);
    plot(B(:, r), 'LineWidth', 1.5);
    title(sprintf('Emission %d', r));
    ylabel('Intensity');

    subplot(3, R, 2*R + r);
    plot(C(:, r), 'LineWidth', 1.5);
    title(sprintf('Excitation %d', r));
    ylabel('Intensity');
end
sgtitle('Best CP Decomposition');
