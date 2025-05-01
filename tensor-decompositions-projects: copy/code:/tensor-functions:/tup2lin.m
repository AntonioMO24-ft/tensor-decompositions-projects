function [alpha] = tup2lin(point, dimensions)
    % Function to convert a tuple into its respective linear index
    % -----------------------------
    % Inputs: point - valid tuple as an array. dimensions - dimensions of
    % the tensor as an array.
    % -----------------------------
    % Outputs: alpha - point's respective linear index.
    
    answer = input('Natural or Reverse? (enter n or r): ', 's');
    if answer == 'n'
        s_1 = 1; % first stride is always 1 for the natural ordering
        n = length(dimensions);
        strides = zeros(n,1);
        strides(1) = s_1;
        % next for loop is to calculate the strides of the given tensor. The
        % strides formulas are given also in Definition 2.2. This saves the
        % strides as a list.
        for i = 2:n
            s_next = s_1*dimensions(i-1);
            strides(i) = s_next;
            s_1 = s_next;
        end
    % new lines from the last time: we added an option to compute
    % everything but now in reverse linear indexing if the user wants to:
    elseif answer == 'r'
        s_n = 1;
        n = length(dimensions);
        strides = zeros(n,1);
        strides(n) = s_n;
        for i = n:-1:2
            s_previous = s_n*dimensions(i);
            strides(i-1) = s_previous;
            s_n = s_previous;
        end
    end
    alpha = 1;
    % next for loop calculates alpha
    for i = 1:n
        x = strides(i)*(point(i) - 1);
        alpha = alpha + x;
    end
end