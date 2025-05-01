% Script to plot all the factor matrices of the rank-3 CP decompositions
% that were calculated in the EEM_elbow_method_rank.m script.

load('rank_solutions'); % loading our results from 9.8

ranks = [3, 4, 5]; % ranks that we are going to use
initializations = 5; 
num_components = size(rank_solutions{3, 1}.U{1}, 2); % Number of components in the rank-3 model should be 3

for trial = 1:initializations
    M = rank_solutions{3, trial}; % Extract CP decomposition for rank-3 (trial)
    
    % factor matrices
    A = M.U{1}; % Mode-1 (Samples)
    B = M.U{2}; % Mode-2 (Emission)
    C = M.U{3}; % Mode-3 (Excitation)
    
    % figure for each initialization
    figure;
    sgtitle(sprintf('Best Rank-3 CP Decomposition - Trial %d', trial));

    % plotting all components
    for comp = 1:num_components
        % samples
        subplot(3, num_components, comp);
        bar(A(:, comp)); % Each column of A represents one component
        title('Sample');
        
        % emission
        subplot(3, num_components, num_components + comp);
        plot(B(:, comp), '-o', 'LineWidth', 1.5);
        title('Emission');
        
        % excitation
        subplot(3, num_components, 2*num_components + comp);
        plot(C(:, comp), '-o', 'LineWidth', 1.5);
        title('Excitation');
    end
end





