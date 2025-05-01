function IP = KruskalIP(Xsubs, Xvals, A, B, C)
    % Function that computes ⟨X, [A, B, C]⟩ efficiently
    % --------------------
    % Inputs: Xsubs - N×3 matrix of subscripts of nonzeros in X 
    % (1-based indexing), Xvals - N×1 vector of values at those subscripts
    % A - m×R factor matrix, B - n×R factor matrix, C - p×R factor matrix.
    % --------------------
    % Outputs: IP - scalar inner product ⟨X, [A, B, C]⟩

    N = length(Xvals);           % number of nonzeros
    %R = size(A, 2);              % rank of Kruskal tensor
    IP = 0;                      % initialize inner product

    for t = 1:N
        i = Xsubs(t, 1);         % mode-1 index
        j = Xsubs(t, 2);         % mode-2 index
        k = Xsubs(t, 3);         % mode-3 index
        x = Xvals(t);            % value of X at (i,j,k)

        % Get the i-th, j-th, and k-th rows from factor matrices
        a = A(i, :);             % 1×R
        b = B(j, :);             % 1×R
        c = C(k, :);             % 1×R

        % Compute contribution: sum(a .* b .* c)
        contrib = sum(a .* b .* c);

        % Accumulate weighted contribution
        IP = IP + x * contrib;
    end
end
