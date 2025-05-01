% Script that uses an existing nonnegative CP-Decomposition applied to the
% Chicago crime data tensor. The rank used is 7.

% Nonnegative CP decomposition
rank = 7;
M = cp_nmu(X, rank);

% Extract factor matrices
A = M.U{1}; % Days (365 × 7)
B = M.U{2}; % Hours (24 × 7)
C = M.U{3}; % Communities (77 × 7)
D = M.U{4}; % Crime Types (12 × 7)

% Plot settings
component_groups = {1:3, 4:6, 7};
for fig_idx = 1:length(component_groups)
    figure;
    components = component_groups{fig_idx};
    num_components = length(components);
    
    for i = 1:num_components
        r = components(i);
        % Row index in subplot
        subplot(num_components, 4, (i-1)*4 + 1);
        plot(A(:, r), 'LineWidth', 1.5);
        title('Days');
        
        subplot(num_components, 4, (i-1)*4 + 2);
        bar(B(:, r));
        title('Hours');
       
        subplot(num_components, 4, (i-1)*4 + 3);
        bar(C(:, r));
        title('Communities');
        
        subplot(num_components, 4, (i-1)*4 + 4);
        bar(D(:, r));
        title('Types');
    end
    
    sgtitle(['CP-NMU Decomposition - Components ', num2str(components)]);
end

