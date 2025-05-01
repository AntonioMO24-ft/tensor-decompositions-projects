function [f, g] = cp_fg_wrapper(v, X, chi, r)
    [f, g] = CP_FG(X, chi, v, r);
end
