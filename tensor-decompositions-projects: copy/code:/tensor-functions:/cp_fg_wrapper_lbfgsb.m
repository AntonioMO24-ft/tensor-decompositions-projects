function [f, g] = cp_fg_wrapper_lbfgsb(v, X, chi, r, nu)
    % Wrapper function so that we can use v only as an input later and also
    % to apply regularizations to f and g.
    [f_cp, g_cp] = CP_FG(X, chi, v, r);

    % Regularizations

    f = f_cp + nu * (v' * v);
    g = g_cp + 2 * nu * v;
end
