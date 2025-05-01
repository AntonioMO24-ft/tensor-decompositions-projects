function [tuples] = lin2tup(alpha, dimensions)
    % Function to convert a linear index into its respective tuple.
    % Here we are converting a linear index to its tuple by using the strides
    % and index formula given in Definition 2.2. 
    % Inputs: alpha - a valid linear index. dimensions - an array containing
    % the dimensions of the tensor.
    % -----------------------------
    % Outputs: tuples - alpha's tuple notation. 

    answer = input('Natural or Reverse? (enter n or r): ', 's');
    if answer == 'n'
        s_1 = 1; % first stride is always 1 in natural ordering
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
    % new lines to make so that the user lets the code know if they want
    % the answer in natural or reverse ordering
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
    tuples = zeros(n,1);
    % next for loop calculates the tuple
    for i = 1:n
        tuples(i) = 1 + floor(mod(alpha - 1, dimensions(i)*strides(i))/strides(i));
    end
end

% expansion: print something for non valid alpha, i.e. bigger than product
% of the dimensions. return the answer as a text as well