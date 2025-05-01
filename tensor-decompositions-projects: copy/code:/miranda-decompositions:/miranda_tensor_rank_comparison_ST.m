% Script to compute the Tucker decomposition of the Miranda tensor with
% ST_HOSVD instead of HOSVD. The first part is for a (600,100,100) rank and
% the second part is for a (300,50,50) rank.
rank_values = [600, 100, 100];
tic;
[G, U, V, W, ERR] = ST_HOSVD(density, rank_values);
elapsed_time = toc;
fprintf('HOSVD took %.4f seconds.\n', elapsed_time);
fprintf('Approximation error: %.6f\n', ERR);

% Output:
% ST-HOSVD took 195.9249 seconds.
% Approximation error: 0.445869
% Took less time than HOSVD and around the same approximation error

rank_values = [300, 50, 50];
tic;
[G, U, V, W, ERR] = ST_HOSVD(density, rank_values);
elapsed_time = toc;
fprintf('ST-HOSVD took %.4f seconds.\n', elapsed_time);
fprintf('Approximation error: %.6f\n', ERR);

% Output:
% ST-HOSVD took 187.3490 seconds.
% Approximation error: 0.463573
% Again took less time than HOSVD and around the same approximation error