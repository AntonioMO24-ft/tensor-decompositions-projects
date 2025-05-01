% Script to find what permutation best matches Atrue to the A in the best 
% rank-3 tensor computed in the EEM_elbow_method_rank.m script. Fortunately
% the tensor has true mixtures of the chemicals that are stored in Atrue
% here. We also find the cosines of the angles between the matched vectors.
% If the best rank-3 approximation was actually good we should see cosines
% of angles close to 1.

% the esimated A is the best rank 3 estimation matrix's first factor matrix
A_est = best_M.U{1};

% this is the given mixture matrix
Atrue = [5.00 0.00 0.00;
         0.00 5.00 0.00;
         0.00 0.00 5.00;
         1.25 5.00 3.75;
         3.75 1.25 5.00;
         5.00 3.75 2.50;
         3.75 3.75 5.00;
         6.25 1.25 1.25;
         1.25 5.00 2.50;
         2.50 6.25 2.50;
         5.00 1.25 3.75;
         1.25 3.75 2.50;
         2.50 3.75 1.25;
         3.75 0.00 2.50;
         2.50 0.00 3.75;
         5.00 0.00 1.25;
         3.75 0.00 3.75;
         3.75 0.00 5.00];

% here we are finding the best permutation
perm = perms(1:3);
num_perms = size(perm, 1);
min_diff = inf;
best_perm = [];

for i = 1:num_perms
    A_perm = A_est(:, perm(i, :)); 
    diff = norm(A_perm - Atrue, 'fro'); 
    if diff < min_diff
        min_diff = diff;
        best_perm = perm(i, :); % and storing the best permuatation
    end
end

% now we apply the permutation to the estimated A
A_best = A_est(:, best_perm);

% and computing the cosine of the angles of the matched vectors
cos_angles = zeros(size(Atrue, 1), 1);
for i = 1:size(Atrue, 1)
    cos_angles(i) = dot(Atrue(i, :), A_best(i, :)) / (norm(Atrue(i, :)) * norm(A_best(i, :)));
end

disp(best_perm); % prints the best permutation
disp(cos_angles); % and the cosine angles


% Output:
% Best column permutation:
%      3     2     1
% 
% Cosines of angles between matched vectors:
%     0.9998
%     0.9992
%     0.9992
%     0.9996
%     0.9990
%     0.9967
%     0.9979
%     0.9954
%     0.9993
%     0.9992
%     0.9982
%     0.9981
%     0.9961
%     0.9935
%     0.9988
%     0.9959
%     0.9975
%     0.9997
% 
% so cosine of the angles are really close to 1 which means the angles are
% really close to 0 so they are closley matched


