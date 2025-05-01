% Script to implement HOSVD on the Miranda tensor a rank of 232×43×41 
% timing each call to LLSV and the TTM compression step and comparing their
% timings. We use the function HOSVD3_timed

rank_values = [232,43,41];
[G, U, V, W, ERR, timings] = HOSVD3_timed(density, rank_values);
% Output: 
% LLSV timings: Mode-1: 145.1797, Mode-2: 30.0154, Mode-3: 31.0044
% TTM timings: Mode-1: 27.0818, Mode-2: 7.4027, Mode-3: 10.7623

% here the mode-1 significantly takes longer since the dimension has
% greater size. The mode-2 and 3 are comparable since they are similar in
% size and only differe by a second on the LLSV and 3 seconds on the TTM
% function.