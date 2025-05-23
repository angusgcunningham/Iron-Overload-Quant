%clc
%clear all
%close all

T2 = [5, 10, 15, 20];       % T2* relaxation times
TE = [1:1.375:16.5]';  % Time to eecho
%TE = 1.375;  % try this
s0 = [155, 255, 355, 455];  % Initial signal inensity values
sigma2 = 0;                 % noise variance
Nrow = 32;                  % Phantom image size
Ncol = 32;                  % Phantom image size

%% Create phantom
for uu = 1:length(T2)
    Phantom_WO_NoiseTemp{uu} = createPhantoms('exp', TE, T2(uu), s0(uu), Nrow, Ncol);
end
Phantom_WO_Noise = [Phantom_WO_NoiseTemp{1}, Phantom_WO_NoiseTemp{2}; Phantom_WO_NoiseTemp{3}, Phantom_WO_NoiseTemp{4}];
[Nrow_, Ncol_, bands] = size(Phantom_WO_Noise);

%% Add noise if required
Y = Phantom_WO_Noise + sqrt(sigma2) * randn(Nrow_, Ncol_, bands);
yReshaped = reshape(Y, Nrow_*Ncol_, bands)';

%% Run T2* estimation algorithm
lambdaA = 1e-5;
lambdaR = 1e-5;
[a, r] = relaxationEst(yReshaped, TE, Nrow_, Ncol_, lambdaA, lambdaR);

%% Plot stuff - MODIFY AS REQUIRED
figure(1)
subplot(2, 4, 1);
% Ground truth s0
S0_Image = [s0(1)*ones(Nrow, Ncol), s0(2)*ones(Nrow, Ncol); s0(3)*ones(Nrow, Ncol), s0(4)*ones(Nrow, Ncol)];
imagesc(S0_Image)
axis image
axis off
caxis([min(min(Y(:, :,1))) max([max(a(:)) s0(uu)])])
c = colorbar;
set(c, 'FontSize', 26)
title('a_0')

subplot(2, 4, 2);
imagesc(Y(:, :, 1))
axis image
axis off
caxis([min(min(Y(:, :,1))) max([max(a(:)) s0(uu)])])
c = colorbar;
set(c, 'FontSize', 26)
title('a_1')

subplot(2, 4, 3)
imagesc(reshape(a(:, :, 1), Nrow_, Ncol_))
axis image
axis off
caxis([min(min(Y(:, :,1))) max([max(a(:)) s0(uu)])])
c = colorbar;
set(c, 'FontSize', 26)
title('Estimated a_0')

a_reshaped = reshape(a, Nrow_, Ncol_);
a_mean = [mean(mean(a_reshaped(1:32, 1:32)))*ones(Nrow, Ncol), mean(mean(a_reshaped(1:32, 33:end)))*ones(Nrow, Ncol); mean(mean(a_reshaped(33:end, 1:32)))*ones(Nrow, Ncol), mean(mean(a_reshaped(33:end, 33:end)))*ones(Nrow, Ncol)];

subplot(2, 4, 4)
imagesc(a_mean)
axis image
axis off
caxis([min(min(Y(:, :,1))) max([max(a_mean(:)) s0(uu)])])
c = colorbar;
set(c, 'FontSize', 26)
title('Mean of Estimated a_0')


figure(2)
subplot(2, 4, 1)
imagesc(reshape(1./r(:, :, 1), Nrow_, Ncol_))
axis image
axis off
caxis([0 max([max(1./r), max(T2)])])
c = colorbar;
set(c, 'FontSize', 26)
colormap hsv
c.TickLabels = [0, 5, 10, 15, 20];
title('Estimated T2*')

r_reshaped = reshape(1./r, Nrow_, Ncol_);
r_mean = [mean(mean(r_reshaped(1:32, 1:32)))*ones(Nrow, Ncol), mean(mean(r_reshaped(1:32, 33:end)))*ones(Nrow, Ncol); mean(mean(r_reshaped(33:end, 1:32)))*ones(Nrow, Ncol), mean(mean(r_reshaped(33:end, 33:end)))*ones(Nrow, Ncol)];

subplot(2, 4, 2)
imagesc(r_mean)
axis image
axis off
caxis([0 max([max(r_mean(:)), max(T2)])])
c = colorbar;
set(c, 'FontSize', 26)
colormap hsv
c.TickLabels = [0, 5, 10, 15, 20];
title('Mean Estimated T2*')


