% Script for running CP-OPT using L-BFGS-B on the EEM tensor 5 times, saving the 
% best solution, according to the relative error and plotting it. Here we
% use l = 0 (enforcing nonnegativity on v) with a regularization of nu = 0.
% This solution should be the same as CP-ALS since they are solving the
% same optimization problem.

% Best error: 0.0331

clear
load('EEM18.mat')

% getting the dimensions of the tensor
dimensions = size(X);
m = dimensions(1);
n = dimensions(2);
p = dimensions(3);


r = 3; % CP rank
num_trials = 5;

errors = zeros(num_trials, 1);       % Store relative errors
solutions = cell(num_trials, 1);     % Store CP decompositions
chi = norm(X)^2;                 % Frobenius norm squared

nu = 0; % regularization parameter
fun = @(v) cp_fg_wrapper_lbfgsb(v, X, chi, r, nu);

%% Main Loop
for i = 1:num_trials
    rng(i); 

    % Initial guess
    v0 = randn(r * (m + n + p), 1);
    
    nVars = length(v0);
    lb = zeros(nVars, 1);   % enforcing v ≥ 0
    ub = inf(nVars, 1);

    % Options
    options = struct();
    options.x0 = v0;
    options.m = 5;                         
    options.maxIts = 1000;                
    options.maxTotalIts = 10000;          
    options.pgtol = 1e-5;                 
    options.factr = 1e-9 / eps;           
    options.printEvery = 10;              
    options.verbose = 1;                  

    [v_opt, fval, info] = lbfgsb(fun, lb, ub, options);


    A = reshape(v_opt(1:m*r), m, r);
    B = reshape(v_opt(m*r+1:m*r+n*r), n, r);
    C = reshape(v_opt(m*r+n*r+1:end), p, r);

    % Create ktensor from factors and weights
    M = ktensor({A, B, C});

    M = arrange(M);

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
figure('Name', 'Best CP-OPT Decomposition (l = 0, nu = 0)');
for r = 1:R
    subplot(3, R, r);
    bar(A(:, r), "red");
    title(sprintf('Sample %d', r));
    ylabel('Weight');

    subplot(3, R, R + r);
    plot(B(:, r), "red", 'LineWidth', 1.5);
    title(sprintf('Emission %d', r));
    ylabel('Intensity');

    subplot(3, R, 2*R + r);
    plot(C(:, r), "red", 'LineWidth', 1.5);
    title(sprintf('Excitation %d', r));
    ylabel('Intensity');
end
sgtitle('Best CP-OPT Decomposition (l = 0, nu = 0)');
