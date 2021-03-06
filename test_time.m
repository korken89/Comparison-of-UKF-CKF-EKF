%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare execution time of the SR-EKF, SR-UKF and SR-CKF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

fprintf('\n')
fprintf('---------------------------------------------------------------')
fprintf('\n')
fprintf('-   Compare execution time of the SR-UKF, SR-CKF and SR-EKF   -')
fprintf('\n')
fprintf('-        By Emil Fresk, Lulea University of Technology        -')
fprintf('\n')
fprintf('---------------------------------------------------------------')
fprintf('\n\n')

% Number of iterations
iter = 100000;

% State size
state_s = 15;

% Measuremenat size
meas_s_best = 6;
meas_s_worst = 9;

% % % % % % % % % % % % % % % % % % % % % %
% Initialize matricies
% % % % % % % % % % % % % % % % % % % % % %

x = zeros(state_s, 1);
u = 0;
z = zeros(meas_s_best, 1);
Sp = eye(state_s);
Sq = eye(state_s);
Sr = eye(meas_s_best);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SR-UKF test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running SR-UKF best case test...')


tic
for i=1:iter
  [ x, Sp ] = srukf_innovate( x, u, z, Sp, Sq, Sr, @f_test, @h_best_test);
end
ukf_time_best = toc;

fprintf('   Done!\n')

% % % % % % % % % % % % % % % % % % % % % %
% Initialize matricies
% % % % % % % % % % % % % % % % % % % % % %

x = zeros(state_s, 1);
u = 0;
z = zeros(meas_s_worst, 1);
Sp = eye(state_s);
Sq = eye(state_s);
Sr = eye(meas_s_worst);

fprintf('Running SR-UKF worst case test...')

tic
for i=1:iter
  [ x, Sp ] = srukf_innovate( x, u, z, Sp, Sq, Sr, @f_test, @h_worst_test);
end
ukf_time_worst = toc;

fprintf('  Done!\n')

% % % % % % % % % % % % % % % % % % % % % %
% Initialize matricies
% % % % % % % % % % % % % % % % % % % % % %

x = zeros(state_s, 1);
u = 0;
z = zeros(meas_s_best, 1);
Sp = eye(state_s);
Sq = eye(state_s);
Sr = eye(meas_s_best);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SR-UKF test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running SR-CKF best case test...')

tic
for i=1:iter
  [ x, Sp ] = srckf_innovate( x, u, z, Sp, Sq, Sr, @f_test, @h_best_test);
end
ckf_time_best = toc;

fprintf('   Done!\n')

% % % % % % % % % % % % % % % % % % % % % %
% Initialize matricies
% % % % % % % % % % % % % % % % % % % % % %

x = zeros(state_s, 1);
u = 0;
z = zeros(meas_s_worst, 1);
Sp = eye(state_s);
Sq = eye(state_s);
Sr = eye(meas_s_worst);

fprintf('Running SR-CKF worst case test...')

tic
for i=1:iter
  [ x, Sp ] = srckf_innovate( x, u, z, Sp, Sq, Sr, @f_test, @h_worst_test);
end
ckf_time_worst = toc;

fprintf('  Done!\n')

% % % % % % % % % % % % % % % % % % % % % %
% Initialize matricies
% % % % % % % % % % % % % % % % % % % % % %

x = zeros(state_s, 1);
u = 0;
z = zeros(meas_s_best, 1);
Sp = eye(state_s);
Sq = eye(state_s);
Sr = eye(meas_s_best);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run SR-UKF test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running SR-EKF best case test...')

tic
for i=1:iter
  [ x, Sp ] = srekf_innovate_best( x, u, z, Sp, Sq, Sr, @f_test, @h_best_test);
end
ekf_time_best = toc;

fprintf('   Done!\n')

% % % % % % % % % % % % % % % % % % % % % %
% Initialize matricies
% % % % % % % % % % % % % % % % % % % % % %

x = zeros(state_s, 1);
u = 0;
z = zeros(meas_s_worst, 1);
Sp = eye(state_s);
Sq = eye(state_s);
Sr = eye(meas_s_worst);

fprintf('Running SR-EKF worst case test...')

tic
for i=1:iter
  [ x, Sp ] = srekf_innovate_worst( x, u, z, Sp, Sq, Sr, @f_test, @h_worst_test);
end
ekf_time_worst = toc;

fprintf('  Done!\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['Iterations:   ', num2str(iter)])
fprintf('\n\n')
fprintf('           Best (s)\t\tWorst (s)\t\tRelative best\t\tRelative worst')
fprintf('\n')
fprintf('           --------\t\t---------\t\t-------------\t\t--------------')
fprintf('\n')
fprintf(['SR-UKF:    ', num2str(ukf_time_best, 3), '  \t\t', num2str(ukf_time_worst, 3)])
fprintf(['  \t\t\t', num2str(ukf_time_best/ekf_time_best, 3), '\t\t\t\t', num2str(ukf_time_worst/ekf_time_worst, 3)])
fprintf('\n')
fprintf(['SR-CKF:    ', num2str(ckf_time_best, 3), '  \t\t', num2str(ckf_time_worst, 3)])
fprintf(['  \t\t\t', num2str(ckf_time_best/ekf_time_best, 3), '\t\t\t\t', num2str(ckf_time_worst/ekf_time_worst, 3)])
fprintf('\n')
fprintf(['SR-EKF:    ', num2str(ekf_time_best, 3), '  \t\t', num2str(ekf_time_worst, 3)])
fprintf(['  \t\t\t', num2str(ekf_time_best/ekf_time_best, 3), '     \t\t\t\t', num2str(ekf_time_worst/ekf_time_worst, 3)])
fprintf('\n\n')