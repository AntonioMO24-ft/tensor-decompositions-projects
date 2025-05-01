% Scripts that compares the relative error of finding the rank-3 cp
% decomposition using als, opt, and nls methods. The resulting figures
% using 5 and 20 random initializations are saved as EEM_cp_comparison5.fig
% and EEM_cp_comparison20.fig, respectively in the figures folder.

clear
load('EEM18.mat')

r = 3; % CP rank
num_trials = 20;

% getting the dimensions of the tensor
dimensions = size(X);
m = dimensions(1);
n = dimensions(2);
p = dimensions(3);

chi = norm(X)^2;

opts = struct('max_iter', 500, 'tol', 1e-10);

nu = 0; % regularization parameter
fun = @(v) cp_fg_wrapper_lbfgsb(v, X, chi, r, nu);

% preallocating the error arrays for each cp alternative
errors_als = zeros(num_trials, 1);       % Store relative errors
errors_opt = zeros(num_trials, 1);
errors_nls = zeros(num_trials, 1);

for i = 1:num_trials
    rng(i); 

    %% Running CP-ALS
    fprintf('Running CP-ALS...\n');
    [A_als, B_als, C_als, lambda] = CP_ALS3(X, r, 50, 1e-4);

    % Create ktensor from factors and weights
    M_als = ktensor(lambda(:), {A_als, B_als, C_als});

    M_als = arrange(M_als);
    % Fix the signs. The key to having nonnegative are these two lines
    
    M_als = fixsigns(M_als);

    % Compute approximation error
    X_approx_als = full(M_als);  % reconstruct full tensor
    errors_als(i) = norm(X(:) - X_approx_als(:)) / norm(X(:));
    
    %% Running CP-OPT
    fprintf('Running CP-OPT...\n');
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


    A_opt = reshape(v_opt(1:m*r), m, r);
    B_opt = reshape(v_opt(m*r+1:m*r+n*r), n, r);
    C_opt = reshape(v_opt(m*r+n*r+1:end), p, r);

    % Create ktensor from factors and weights
    M_opt = ktensor({A_opt, B_opt, C_opt});

    M_opt = arrange(M_opt);

    % Compute approximation error
    X_approx_opt = full(M_opt);  % Converts to full tensor
    errors_opt(i) = norm(X(:) - X_approx_opt(:)) / norm(X(:));

    %% Running CP-NLS  
    fprintf('Running CP-NLS...\n');
    A0 = randn(m, r);
    B0 = randn(n, r);
    C0 = randn(p, r);
    v_init = [A0(:); B0(:); C0(:)];

    % Run CP-NLS
    [v_final, output] = CP_NLS(v_init, X, chi, dimensions, r, opts);

    A_nls = reshape(v_final(1:m*r), m, r);
    B_nls = reshape(v_final(m*r+1:m*r+n*r), n, r);
    C_nls = reshape(v_final(m*r+n*r+1:end), p, r);

    % Create ktensor from factors and weights
    M_nls = ktensor({A_nls, B_nls, C_nls});

    M_nls = arrange(M_nls);

    M_nls = fixsigns(M_nls);

    % Compute approximation error
    X_approx_nls = full(M_nls);  % Converts to full tensor
    errors_nls(i) = norm(X(:) - X_approx_nls(:)) / norm(X(:));
end

%% Plot comparison
figure;
hold on;
plot(0:length(errors_als)-1, errors_als, '-o', 'DisplayName', 'CP-ALS');
plot(0:length(errors_nls)-1, errors_nls, '-x', 'DisplayName', 'CP-NLS');
plot(0:length(errors_opt)-1, errors_opt, '-s', 'DisplayName', 'CP-OPT');
xlabel('Iteration');
ylabel('Relative Error');
legend('show');
title('Convergence Comparison of CP-ALS, CP-NLS, and CP-OPT');
grid on;
