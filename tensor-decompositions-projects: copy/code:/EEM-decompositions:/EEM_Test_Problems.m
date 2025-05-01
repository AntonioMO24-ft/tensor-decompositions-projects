% Script to create test problems with known rank, with and without noise, to test 
% the method. Using both 3- and 4-way tensors, of sizes 
% 100 × 80 × 60 and r = 10 and 50 × 40 × 30 × 20 and r = 8.using random 
% factor matrices, e.g., using the Tensor Toolbox command 
% A=matrandnorm(n,r) as well as matrices that have congruent factors via 
% A=matrandcong(n,r,gamma), which ensures that (a^T_i a_j)/(||a_i||2||a_j||_2) = γ 
% for all pairs of columns of A. For γ = 0, all columns are orthogonal. 
% For γ = 1, all columns are identical. The problems generally become more 
% difficult for γ ≥ 0.5.

clear
load('EEM18.mat')

dims = [100, 80, 60];
rank = 10;

% random factor matrices
A = matrandnorm(dims(1), rank);
B = matrandnorm(dims(2), rank);
C = matrandnorm(dims(3), rank);

% Constructing the tensor via ktensor
lambda_true = ones(rank, 1);
X_true = ktensor(lambda_true, {A, B, C});
X = full(X_true);

% Adding noise
% noise_level = 0.1;  % 10% noise
% X = X + noise_level * tensor(randn(size(X)));

% Running CP-ALSd
tol = 1e-10;
max_iters = 50;
[factors_est, lambda_est] = CP_ALSd(X, rank, tol, max_iters);


lambda_est = lambda_est(:);

% Reconstructing the estimated tensor from the CPD
X_est = full(ktensor(lambda_est, factors_est));

% Compute relative error
rel_error = norm(X(:) - X_est(:)) / norm(X(:));
fprintf('Relative error: %.4f\n', rel_error);
