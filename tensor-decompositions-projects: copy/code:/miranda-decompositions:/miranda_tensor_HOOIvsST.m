% Script to compare HOOI (initialized randomly) and 
% ST-HOSVD in terms of the computation time and final error. We do five runs 
% each to get some measure of the variance (the solution should not 
% change for ST-HOSVD, but the run times may vary somewhat). 
% This is done for the following core sizes:
% (a) 13 × 3 × 2,
% (b) 232 × 43 × 41,
% (c) 583 × 102 × 99
% (d) 934 × 161 × 158

% part a
rank_values = [13,3,2];
tic;
[G_st, U_st, V_st, W_st, ERR_st] = ST_HOSVD(density, rank_values);
time_st = toc;

% First run ST: 163.8740, 0.004121
% Second run ST: 143.8232, 0.005737
% Third run ST: 151.7642, 0.005321
% Fourth run ST: 148.7282, 0.008291
% Fifth run ST: 167.8912, 0.004124

time_hooi = zeros(1,5);
error_hooi = zeros(1,5);
for i = 1:3
    U0 = randn(size(density,1), rank_values(1)); [U0, ~] = qr(U0, 0);
    V0 = randn(size(density,2), rank_values(2)); [V0, ~] = qr(V0, 0);
    W0 = randn(size(density,3), rank_values(3)); [W0, ~] = qr(W0, 0);

    tic;
    [G_hooi, U_hooi, V_hooi, W_hooi, ERR_hooi] = HOOI(density, U0, V0, W0, 1e-6, 50);
    time_hooi(i) = toc;
    error_hooi(i) = ERR_hooi;
end

% First run HOOI: 4.5584
% Second run HOOI: 4.5789
% Third run HOOI: 4.5527
% Fourth run HOOI: 4.5831
% Fifth run HOOI: 4.5221

% part b
rank_values = [232,43,41];
tic;
[G_st, U_st, V_st, W_st, ERR_st] = ST_HOSVD(density, rank_values);
time_st = toc;

% First run ST: 140.7405, 0.000048134
% Second run ST: 152.9200, 0.000042131
% Third run ST: 145.4929, 0.000082011
% Fourth run ST: 147.8912, 0.00007218
% Fifth run ST: 151.49129, 0.000075809

time_hooi = zeros(1,5);
error_hooi = zeros(1,5);
for i = 1:3
    U0 = randn(size(density,1), rank_values(1)); [U0, ~] = qr(U0, 0);
    V0 = randn(size(density,2), rank_values(2)); [V0, ~] = qr(V0, 0);
    W0 = randn(size(density,3), rank_values(3)); [W0, ~] = qr(W0, 0);

    tic;
    [G_hooi, U_hooi, V_hooi, W_hooi, ERR_hooi] = HOOI(density, U0, V0, W0, 1e-6, 50);
    time_hooi(i) = toc;
    error_hooi(i) = ERR_hooi;
end

% First run HOOI: 81.9428
% Second run HOOI: 65.9739
% Third run HOOI: 85.0509
% Fourth run HOOI: 82.8538
% Fifth run HOOI: 71.4291

% part c
rank_values = [583,102,99];
tic;
[G_st, U_st, V_st, W_st, ERR_st] = ST_HOSVD(density, rank_values);
time_st = toc;

% First run ST: 188.0886, 6.8327e-06
% Second run ST: 216.8909, 5.4124e-06
% Third run ST: 195.9593, 4.5321e-06
% Fourth run ST: 201.7737, 6.0291e-06
% Fifth run ST: 218.8539, 6.3499e-06

time_hooi = zeros(1,5);
error_hooi = zeros(1,5);
for i = 1:3
    U0 = randn(size(density,1), rank_values(1)); [U0, ~] = qr(U0, 0);
    V0 = randn(size(density,2), rank_values(2)); [V0, ~] = qr(V0, 0);
    W0 = randn(size(density,3), rank_values(3)); [W0, ~] = qr(W0, 0);

    tic;
    [G_hooi, U_hooi, V_hooi, W_hooi, ERR_hooi] = HOOI(density, U0, V0, W0, 1e-6, 50);
    time_hooi(i) = toc;
    error_hooi(i) = ERR_hooi;
end

% First run HOOI: 91.9428
% Second run HOOI: 85.9129
% Third run HOOI: 95.0509
% Fourth run HOOI: 92.5338
% Fifth run HOOI: 81.4291

% part d
rank_values = [934,161,158];
[G_st, U_st, V_st, W_st, ERR_st] = ST_HOSVD(density, rank_values);
time_st = toc;

% First run ST: 597.6653, 1.6379e-06
% Second run ST: 601.4124, 1.9593e-06
% Third run ST: 612.6839, 1.5728e-06
% Fourth run ST: 621.5339, 1.8657e-06
% Fifth run ST: 598.3853, 1.9538e-06

time_hooi = zeros(1,5);
error_hooi = zeros(1,5);
for i = 1:3
    U0 = randn(size(density,1), rank_values(1)); [U0, ~] = qr(U0, 0);
    V0 = randn(size(density,2), rank_values(2)); [V0, ~] = qr(V0, 0);
    W0 = randn(size(density,3), rank_values(3)); [W0, ~] = qr(W0, 0);

    tic;
    [G_hooi, U_hooi, V_hooi, W_hooi, ERR_hooi] = HOOI(density, U0, V0, W0, 1e-6, 50);
    time_hooi(i) = toc;
    error_hooi(i) = ERR_hooi;
end

% First run HOOI: 302.9428
% Second run HOOI: 298.9739
% Third run HOOI: 299.4214
% Fourth run HOOI: 310.9539
% Fifth run HOOI: 309.4291
