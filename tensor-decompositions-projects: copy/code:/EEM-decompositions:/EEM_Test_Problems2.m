% Script to create test problems with known rank, with and without noise, to test the 
% convergence of CP-NLS. Generates three random 250 × 10 factor matrices 
% with specified congruence γ = {0.1, 0.8}, e.g., using the Tensor Toolbox 
% command A=matrandcong(n,r,gamma). If X is the tensor with exact rank 10, 
% compute Z = X + η∥X∥∥N∥· N, where η is the fraction of noise to add and N 
% is a randomly generated tensor. 

clear; clc;

% Settings
gammas = [0.1, 0.8];   % congruence levels
etas = [0, 1e-4, 0.1]; % nise levels
num_trials = 5;   

n = 250;   % size of each mode
r = 10;    % rank
dims = [n, n, n];

opts = struct('max_iter', 500, 'tol', 1e-10); % options for NLS

% Plot
figure;
tiledlayout(length(gammas), length(etas), 'TileSpacing', 'Compact');

trial_colors = lines(num_trials); % colors for trials

for gi = 1:length(gammas)
    gamma = gammas(gi);

    % ground truth factor matrices
    Atrue = matrandcong(n, r, gamma);
    Btrue = matrandcong(n, r, gamma);
    Ctrue = matrandcong(n, r, gamma);

    % reconstructing to a full tensor
    X = ktensor({Atrue, Btrue, Ctrue});
    X = full(X);

    for ei = 1:length(etas)
        eta = etas(ei);

        % adding noise
        N = tensor(randn(size(X)));
        noise_scale = eta * (norm(X) / norm(N));
        Z = X + noise_scale * N;

        % computing the norm
        chi = norm(Z)^2;

        nexttile;
        hold on;
        title(sprintf('\\gamma=%.1f, \\eta=%.4f', gamma, eta));
        xlabel('Iteration');
        ylabel('Relative Error');
        set(gca, 'YScale', 'log');

        % Running multiple trials
        for trial = 1:num_trials
            fprintf('Running trial %d for gamma=%.1f, eta=%.4e\n', trial, gamma, eta);

            % Random initialization
            A0 = randn(n, r);
            B0 = randn(n, r);
            C0 = randn(n, r);
            v_init = [A0(:); B0(:); C0(:)];

            % Running CP-NLS
            [v_final, output] = CP_NLS(v_init, Z, chi, dims, r, opts);

            % Plotting convergence curve
            plot(0:length(output.relerrs)-1, output.relerrs, '-', 'Color', trial_colors(trial,:), 'LineWidth', 1.2);
        end

        hold off;
    end
end

sgtitle('CP-NLS Convergence for Different Noise and Congruence Settings');
