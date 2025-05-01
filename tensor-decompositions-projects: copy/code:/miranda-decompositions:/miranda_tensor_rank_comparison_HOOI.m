% Script to calculate the Tucker Decomposition of the Miranda Tensor for a
% (600,100,100) rank and a (60,10,10) rank tensor and comparing them. Here
% we use the HOOI function instead of calculating the HOSVD of the tensor.

tau = 1e-6;
MAXITERS = 100;

% Part a
ranks = [600, 100, 100];
% Initialize factor matrices (orthonormal)
[U, ~] = qr(randn(size(density,1), ranks(1)), 0);
[V, ~] = qr(randn(size(density,2), ranks(2)), 0);
[W, ~] = qr(randn(size(density,3), ranks(3)), 0);
% Run HOOI
tic;
[G, U, V, W, ERR, num_iters] = HOOI(density, U, V, W, tau, MAXITERS);
time_elapsed = toc;
fprintf('Number of iterations: %d\n', num_iters);
fprintf('Time taken: %.4f seconds\n', time_elapsed);
fprintf('Final error: %.6f\n', ERR);

% Output:
% Number of iterations: 6
% Time taken: 114.2843 seconds
% Final error: 0.023130612

% Part b
ranks = [60, 10, 10];
% Initialize factor matrices (orthonormal)
[U, ~] = qr(randn(size(density,1), ranks(1)), 0);
[V, ~] = qr(randn(size(density,2), ranks(2)), 0);
[W, ~] = qr(randn(size(density,3), ranks(3)), 0);
% Run HOOI
tic;
[G, U, V, W, ERR, num_iters] = HOOI(density, U, V, W, tau, MAXITERS);
time_elapsed = toc;
fprintf('Number of iterations: %d\n', num_iters);
fprintf('Time taken: %.4f seconds\n', time_elapsed);
fprintf('Final error: %.6f\n', ERR);

% Output:
% Number of iterations: 6
% Time taken: 19.0580 seconds
% Final error: 0.023130612

% HOOI took significantly less time and it is more accurate than HOSVD and
% ST_HOSVD, as expected



