% Script that uses an existing code to compute the nonnegative CP of the
% Monkey BMI data tensor for ranks 5 through 15 using at least three
% starting points per run. Also plots the best relative error versus its
% respective rank.

% Parameters
ranks = 5:15; % we are going to do CPD on ranks 5 to 15
initializations = 3;
errors = zeros(11, 1);  % storing the errors to plot later

% Doing CPD on each of the ranks
for i = 1:11
    rank = ranks(i);
    best_error = inf;
    
    % Doing CP for each initialization per rank
    for j = 1:initializations
        % here is where we run the cp algorithm
        P = cp_nmu(X, rank); % nmu is a cpd with nonnegative constraints on the factor matrices
        recon = full(P);  % full() reconstructs P into a full tensor
        error = norm(X(:) - recon(:)) / norm(X(:)); % relative error
        best_error = min(best_error, error); % storing the best error 
        % of all initializations per rank
    end
    errors(i) = best_error; % saving it
end
% Plotting the results
figure;
plot(ranks, errors, '-o', 'LineWidth', 2);
xlabel('CP Rank');
ylabel('Best Relative Error');
title('Best Relative Error vs. CP Rank');
grid on;
