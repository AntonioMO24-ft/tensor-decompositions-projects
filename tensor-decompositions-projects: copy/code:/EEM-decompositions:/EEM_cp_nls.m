% Script to get the rank-3 CP decomposition of the EEM tensor. This script
% uses CP-NLS (nonlinear least squares) to compute the decomposition.
% Best error: 0.0329

clear
load('EEM18.mat')

% getting the dimensions of the tensor
dimensions = size(X);
m = dimensions(1);
n = dimensions(2);
p = dimensions(3);


r = 3; % CP rank
num_trials = 5;
opts = struct('max_iter', 500, 'tol', 1e-10);

errors = zeros(num_trials, 1);       % Store relative errors
solutions = cell(num_trials, 1);     % Store CP decompositions
chi = norm(X)^2;                 % Frobenius norm squared

%% Main Loop
for i = 1:num_trials
    rng(i); 

    A0 = randn(m, r);
    B0 = randn(n, r);
    C0 = randn(p, r);
    v_init = [A0(:); B0(:); C0(:)];

    % Run CP-NLS
    [v_final, output] = CP_NLS(v_init, X, chi, dimensions, r, opts);

    A = reshape(v_final(1:m*r), m, r);
    B = reshape(v_final(m*r+1:m*r+n*r), n, r);
    C = reshape(v_final(m*r+n*r+1:end), p, r);

    % Create ktensor from factors and weights
    M = ktensor({A, B, C});

    M = arrange(M);

    M = fixsigns(M);

    % Compute approximation error
    X_approx = full(M);  % Converts to full tensor
    errors(i) = norm(X(:) - X_approx(:)) / norm(X(:));

    % Store CP model
    solutions{i} = M;
end

% Same code as ALS
[best_error, best_idx] = min(errors);
best_model = solutions{best_idx};

% Extract factors
A = best_model.U{1}; B = best_model.U{2}; C = best_model.U{3};
R = size(A, 2);

% Plot 
figure('Name', 'Best CP-NLS Decomposition (l = 0, nu = 0)');
for r = 1:R
    subplot(3, R, r);
    bar(A(:, r), "green");
    title(sprintf('Sample %d', r));
    ylabel('Weight');

    subplot(3, R, R + r);
    plot(B(:, r), "green", 'LineWidth', 1.5);
    title(sprintf('Emission %d', r));
    ylabel('Intensity');

    subplot(3, R, 2*R + r);
    plot(C(:, r), "green", 'LineWidth', 1.5);
    title(sprintf('Excitation %d', r));
    ylabel('Intensity');
end
sgtitle('Best CP-NLS Decomposition');
