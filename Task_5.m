% Task 5a - Report estimated images of both a and T2* = 1/r
figure(3)
subplot(2, 2, 1);
% Ground truth s0
S0_Image = [s0(1)*ones(Nrow, Ncol), s0(2)*ones(Nrow, Ncol); s0(3)*ones(Nrow, Ncol), s0(4)*ones(Nrow, Ncol)];
imagesc(S0_Image)
axis image
axis off
colorbar
title('Ground Truth a (s0)')

subplot(2, 2, 2);
imagesc(a_reshaped)
axis image
axis off
colorbar
title('Estimated a (s0)')

subplot(2, 2, 3);
% Ground truth T2*
T2_Image = [T2(1)*ones(Nrow, Ncol), T2(2)*ones(Nrow, Ncol); T2(3)*ones(Nrow, Ncol), T2(4)*ones(Nrow, Ncol)];
imagesc(T2_Image)
axis image
axis off
colorbar
title('Ground Truth T2* = 1/r')

subplot(2, 2, 4);
imagesc(r_reshaped)
axis image
axis off
colorbar
title('Estimated T2* = 1/r')



% Task 5b - Compute mean of a and T2* = 1/r in each box and compare with ground truth
% Calculate means in each quadrant
a_mean_q1 = mean(mean(a_reshaped(1:Nrow, 1:Ncol)));
a_mean_q2 = mean(mean(a_reshaped(1:Nrow, Ncol+1:end)));
a_mean_q3 = mean(mean(a_reshaped(Nrow+1:end, 1:Ncol)));
a_mean_q4 = mean(mean(a_reshaped(Nrow+1:end, Ncol+1:end)));

r_mean_q1 = mean(mean(r_reshaped(1:Nrow, 1:Ncol)));
r_mean_q2 = mean(mean(r_reshaped(1:Nrow, Ncol+1:end)));
r_mean_q3 = mean(mean(r_reshaped(Nrow+1:end, 1:Ncol)));
r_mean_q4 = mean(mean(r_reshaped(Nrow+1:end, Ncol+1:end)));

% Display estimated means vs ground truth
fprintf('Task 5b - Results without noise (sigma2 = %d):\n', sigma2);
fprintf('Box 1 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(1), T2(1));
fprintf('Box 1 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_mean_q1, r_mean_q1);
fprintf('Box 1 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_mean_q1-s0(1))/s0(1), 100*(r_mean_q1-T2(1))/T2(1));
fprintf('\n');

fprintf('Box 2 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(2), T2(2));
fprintf('Box 2 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_mean_q2, r_mean_q2);
fprintf('Box 2 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_mean_q2-s0(2))/s0(2), 100*(r_mean_q2-T2(2))/T2(2));
fprintf('\n');

fprintf('Box 3 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(3), T2(3));
fprintf('Box 3 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_mean_q3, r_mean_q3);
fprintf('Box 3 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_mean_q3-s0(3))/s0(3), 100*(r_mean_q3-T2(3))/T2(3));
fprintf('\n');

fprintf('Box 4 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(4), T2(4));
fprintf('Box 4 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_mean_q4, r_mean_q4);
fprintf('Box 4 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_mean_q4-s0(4))/s0(4), 100*(r_mean_q4-T2(4))/T2(4));
fprintf('\n\n');


% Task 5c - Add Gaussian noise and repeat
fprintf('Task 5c - Adding Gaussian noise with zero mean and variance of 10\n');
sigma2 = 10;  % Set noise variance to 10

% Create noisy data
Y_noisy = Phantom_WO_Noise + sqrt(sigma2) * randn(Nrow_, Ncol_, bands);
yReshaped_noisy = reshape(Y_noisy, Nrow_*Ncol_, bands)';

% Run T2* estimation algorithm on noisy data
[a_noisy, r_noisy] = relaxationEst(yReshaped_noisy, TE, Nrow_, Ncol_, lambdaA, lambdaR);

% Reshape results for visualization
a_noisy_reshaped = reshape(a_noisy, Nrow_, Ncol_);
r_noisy_reshaped = reshape(1./r_noisy, Nrow_, Ncol_); % Convert r to T2* = 1/r

% Task 5c - Report estimated images with noise
figure(2)
subplot(2, 2, 1);
imagesc(S0_Image)
axis image
axis off
colorbar
title('Ground Truth a (s0)')

subplot(2, 2, 2);
imagesc(a_noisy_reshaped)
axis image
axis off
colorbar
title('Estimated a (s0) with Noise')

subplot(2, 2, 3);
imagesc(T2_Image)
axis image
axis off
colorbar
title('Ground Truth T2* = 1/r')

subplot(2, 2, 4);
imagesc(r_noisy_reshaped)
axis image
axis off
colorbar
title('Estimated T2* = 1/r with Noise')

% Calculate means in each quadrant for noisy case
a_noisy_mean_q1 = mean(mean(a_noisy_reshaped(1:Nrow, 1:Ncol)));
a_noisy_mean_q2 = mean(mean(a_noisy_reshaped(1:Nrow, Ncol+1:end)));
a_noisy_mean_q3 = mean(mean(a_noisy_reshaped(Nrow+1:end, 1:Ncol)));
a_noisy_mean_q4 = mean(mean(a_noisy_reshaped(Nrow+1:end, Ncol+1:end)));

r_noisy_mean_q1 = mean(mean(r_noisy_reshaped(1:Nrow, 1:Ncol)));
r_noisy_mean_q2 = mean(mean(r_noisy_reshaped(1:Nrow, Ncol+1:end)));
r_noisy_mean_q3 = mean(mean(r_noisy_reshaped(Nrow+1:end, 1:Ncol)));
r_noisy_mean_q4 = mean(mean(r_noisy_reshaped(Nrow+1:end, Ncol+1:end)));

% Display estimated means vs ground truth for noisy case
fprintf('Task 5c - Results with noise (sigma2 = %d):\n', sigma2);
fprintf('Box 1 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(1), T2(1));
fprintf('Box 1 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_noisy_mean_q1, r_noisy_mean_q1);
fprintf('Box 1 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_noisy_mean_q1-s0(1))/s0(1), 100*(r_noisy_mean_q1-T2(1))/T2(1));
fprintf('\n');

fprintf('Box 2 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(2), T2(2));
fprintf('Box 2 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_noisy_mean_q2, r_noisy_mean_q2);
fprintf('Box 2 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_noisy_mean_q2-s0(2))/s0(2), 100*(r_noisy_mean_q2-T2(2))/T2(2));
fprintf('\n');

fprintf('Box 3 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(3), T2(3));
fprintf('Box 3 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_noisy_mean_q3, r_noisy_mean_q3);
fprintf('Box 3 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_noisy_mean_q3-s0(3))/s0(3), 100*(r_noisy_mean_q3-T2(3))/T2(3));
fprintf('\n');

fprintf('Box 4 - Ground truth: a = %.2f, T2* = %.2f ms\n', s0(4), T2(4));
fprintf('Box 4 - Estimated:    a = %.2f, T2* = %.2f ms\n', a_noisy_mean_q4, r_noisy_mean_q4);
fprintf('Box 4 - Error:        a = %.2f%%, T2* = %.2f%%\n', 100*(a_noisy_mean_q4-s0(4))/s0(4), 100*(r_noisy_mean_q4-T2(4))/T2(4));

% Comment on results
fprintf('\n\nTask 5 - Results Summary and Comments:\n');
fprintf('1. Without noise, the estimation algorithm performs well, with errors generally below 1%%.\n');
%fprintf('2. With noise (variance = 10), the estimation accuracy decreases, particularly for T2* values.\n');⁠