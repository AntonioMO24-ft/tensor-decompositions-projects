% Script to plot all the factor matrices of the rank-4 and 5 CP 
% decompositions that were calculated in the EEM_elbow_method_rank.m script.

load('rank_solutions'); 
ranks_to_plot = [4, 5]; % Ranks we want to visualize

for rank_idx = 1:length(ranks_to_plot)
    rank = ranks_to_plot(rank_idx);
    
    % finding best trial which is the smallest error
    [~, best_trial] = min(errors(rank, :)); % errors has all errors of cp computed in exercise 9.8
    M = rank_solutions{rank, best_trial}; % best rank -idx approximation from 9.8

    % factor matrices
    A = M.U{1}; % (Samples)
    B = M.U{2}; % (Emission)
    C = M.U{3}; % (Excitation)

    % plotting
    figure;
    sgtitle(sprintf('Best Rank-%d CP Decomposition', rank));

    for comp = 1:rank
        subplot(3, rank, comp);
        bar(A(:, comp)); 
        title('Sample');
    end

    for comp = 1:rank
        subplot(3, rank, rank + comp);
        plot(B(:, comp), 'LineWidth', 1.5);
        title('Emission');
    end

    for comp = 1:rank
        subplot(3, rank, 2*rank + comp);
        plot(C(:, comp), 'LineWidth', 1.5);
        title('Excitation');
    end
end

% best trials: Trial 4 for rank 4 and Trial 1 for rank 5

