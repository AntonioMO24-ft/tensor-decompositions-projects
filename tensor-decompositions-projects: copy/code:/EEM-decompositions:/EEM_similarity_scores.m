% Script that considering the EEM tensor does the following:
% (a) Computes the rank-3 CP with 50 different starting points and
% determines the model that yields the minimal error. 
% (b) Computes the similarity score of the other 49 solutions with the one 
% that yielded the minimal error. 
% (c) Computes the range of similarity scores. 
% (d) Plots the best solution alongside a solution with a high 
% similarity score. 
% (e) Plots the least similar solution alongside the best. 


rank = 3; % we are going to be using rank 3 for the EEM Tensor
num_trials = 50;

errors = zeros(num_trials, 1); % storing the errors
solutions = cell(num_trials, 1); % and each cp decomposition

% part a is a copy of the last exercise of chapter 9.
% (a) Computing 50 CP decompositions
for i = 1:num_trials
    rng(i); 
    [M, ~] = cp_als(X, rank);
    X_approx = full(M);
    error = norm(X(:) - X_approx(:)) / norm(X(:)); % calculating the errors
    errors(i) = error;
    solutions{i} = M;
end
% now we want to store the solution with the minimum error and its index
[best_error, best_idx] = min(errors);

% (b) Compute similarity scores 
% here we will compute similartiy scores using the score function
similarity_scores = zeros(num_trials, 1);
best_model = solutions{best_idx}; % best_model is the best approximation
% which means it has the minimum error

% here we calculate all the similarity scores and store in the same named
% array
for i = 1:num_trials
    if i ~= best_idx
        similarity_scores(i) = score(best_model, solutions{i});
    else
        similarity_scores(i) = 1; % perfect match with itself
    end
end

% (c) Range of similarity scores
% we create the array other_scores so we can do the comparison of finding
% the most similar and least similar. This is why we have to do the
% adjusting in a few lines below
other_scores = similarity_scores([1:best_idx-1, best_idx+1:end]);
fprintf("Similarity score range [%.4f, %.4f]\n", ...
    min(other_scores), max(other_scores));
% here we have to adjust the indices since other scores excludes the best
% trial
[~, most_sim_idx] = max(other_scores);
if most_sim_idx >= best_idx
    most_sim_idx = most_sim_idx + 1;
end
most_similar = solutions{most_sim_idx}; % storing the most similar solution from solutions
[~, least_sim_idx] = min(other_scores);
if least_sim_idx >= best_idx
    least_sim_idx = least_sim_idx + 1;
end
least_similar = solutions{least_sim_idx}; % likewise the least similar solution

% (d) and (e)
% Plots
% the plots are a copy from exercise 9_12. Did not loop because it was
% giving me empty figures so I did it one by one
% first we plot the best model which has index 11
A = best_model.U{1}; B = best_model.U{2}; C = best_model.U{3};
R = size(A, 2);
figure('Name', 'Best CP Decomposition (11)');
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

% here is the plot of the most similar which has index of 20
A = most_similar.U{1}; B = most_similar.U{2}; C = most_similar.U{3};
R = size(A, 2);
figure('Name', 'Most Similar to Best (20) (0.9956)');
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
sgtitle('Most Similar to Best CP Decomposition');

% plot of the least similar with indexf of 7
A = least_similar.U{1}; B = least_similar.U{2}; C = least_similar.U{3};
R = size(A, 2);
figure('Name', 'Least Similar to Best (7) (0.2302)');
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
sgtitle('Least Similar to Best CP Decomposition');
