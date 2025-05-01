% Exercise 6.6a
rank_values = [600, 100, 100];
tic;
[G, U, V, W, ERR] = HOSVD3(density, rank_values);
elapsed_time = toc;
fprintf('HOSVD took %.4f seconds.\n', elapsed_time);
fprintf('Approximation error: %.6f\n', ERR);

% Output:
% HOSVD took 270.5692 seconds.
% Approximation error: 0.309815

% Exercise 6.6b
rank_values = [300, 50, 50];
tic;
[G, U, V, W, ERR] = HOSVD3(density, rank_values);
elapsed_time = toc;
fprintf('HOSVD took %.4f seconds.\n', elapsed_time);
fprintf('Approximation error: %.6f\n', ERR);

% Output:
% HOSVD took 270.7607 seconds.
% Approximation error: 0.970416
