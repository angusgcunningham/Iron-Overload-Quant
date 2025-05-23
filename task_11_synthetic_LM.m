function task_11_synthetic_LM
    % Task 11: Synthetic Phantom LM Fitting with Fixed Color Scales
    % - Amplitude (a) display range: [150, 450] with parula colormap
    % - T2* display range: [5, 20] ms with hsv colormap

    %% 1) Define phantom parameters
    T2    = [5, 10, 15, 20];          % Ground truth T2* (ms)
    s0    = [155, 255, 355, 455];     % Ground truth initial signal
    TE    = (1:1.375:16.5)';          % Echo times (ms)
    Nq    = 32;                       % Quadrant size
    sigma2= 0;                       % Noise variance

    %% 2) Create synthetic phantom (2×2 grid)
    for k = 1:4
        PhantomTemp{k} = createPhantoms('exp', TE, T2(k), s0(k), Nq, Nq);
    end
    Phantom = [PhantomTemp{1}, PhantomTemp{2}; PhantomTemp{3}, PhantomTemp{4}];
    [Nrow, Ncol, nb] = size(Phantom);

    %% 3) Ground truth maps
    S0_gt = [s0(1)*ones(Nq), s0(2)*ones(Nq);
             s0(3)*ones(Nq), s0(4)*ones(Nq)];
    T2_gt = [T2(1)*ones(Nq), T2(2)*ones(Nq);
             T2(3)*ones(Nq), T2(4)*ones(Nq)];

    %% 4) Add noise and reshape for fitting
    Y = Phantom + sqrt(sigma2)*randn(size(Phantom));
    Ymat = reshape(Y, Nrow*Ncol, nb)';  % [nb x (Nrow*Ncol)]

    %% 5) Levenberg-Marquardt fitting pixel-wise
    Npix = Nrow*Ncol;
    a_lm = nan(Npix,1);
    r_lm = nan(Npix,1);
    lb = [0; 0]; ub = [Inf; 0.5];       % r ≤ 0.5 => T2* ≥ 2 ms
    opts = optimset('Display','off');
    model = @(p,t) p(1)*exp(-p(2)*t);
    for pix = 1:Npix
        sig = Ymat(:,pix);
        init = [sig(1), 0.05];
        pp = lsqcurvefit(model, init, TE, sig, lb, ub, opts);
        a_lm(pix) = pp(1);
        r_lm(pix) = pp(2);
    end
    a_lm_img   = reshape(a_lm,   Nrow, Ncol);
    T2_lm_img  = reshape(1./max(r_lm,eps), Nrow, Ncol);

    %% 6) Median filter to reduce speckle
    a_lm_img_filt  = medfilt2(a_lm_img, [3 3]);
    T2_lm_img_filt = medfilt2(T2_lm_img,[3 3]);

    %% 7) Compute quadrant means
    quad_a_lm  = [mean(mean(a_lm_img_filt(1:Nq,   1:Nq))),   mean(mean(a_lm_img_filt(1:Nq,   Nq+1:end)));
                  mean(mean(a_lm_img_filt(Nq+1:end, 1:Nq))),   mean(mean(a_lm_img_filt(Nq+1:end, Nq+1:end)))];
    quad_T2_lm = [mean(mean(T2_lm_img_filt(1:Nq,   1:Nq))),   mean(mean(T2_lm_img_filt(1:Nq,   Nq+1:end)));
                  mean(mean(T2_lm_img_filt(Nq+1:end, 1:Nq))),   mean(mean(T2_lm_img_filt(Nq+1:end, Nq+1:end)))];

    %% 8) Display ground truth vs LM with fixed color scales
    figure('Name','Synthetic Data: GT vs LM','Position',[100 100 1200 600]);

    % Amplitude row
    subplot(2,2,1);
    imagesc(S0_gt, [150 450]);
    colormap(parula);
    colorbar;
    title('Ground Truth a (S0)');
    axis image off;

    subplot(2,2,2);
    imagesc(a_lm_img_filt, [150 450]);
    colormap(parula);
    colorbar;
    title('Estimated a');
    axis image off;

    % T2* row
    subplot(2,2,3);
    imagesc(T2_gt, [5 20]);
    colormap(parula);
    colorbar;
    title('Ground Truth T* = 1/r');
    axis image off;

    subplot(2,2,4);
    imagesc(T2_lm_img_filt, [5 20]);
    colormap(parula);
    colorbar;
    title('Estimated T2*');
    axis image off;

    %% 9) Print quadrant-wise results
    fprintf('\nQuadrant | True a  Est. a  | True T2*  Est. T2*\n');
    fprintf('---------|-----------------|-------------------\n');
    for q = 1:4
        fprintf('   %d     |  %6.1f  %6.1f  |   %6.1f   %6.1f\n', ...
            q, s0(q), quad_a_lm(q), T2(q), quad_T2_lm(q));
    end

    % Error metrics
    a_err  = abs(quad_a_lm  - s0(:)') ./ s0(:)'  * 100;
    T2_err = abs(quad_T2_lm - T2(:)') ./ T2(:)' * 100;
    fprintf('\nMean a error:  %.2f%%\n', mean(a_err));
    fprintf('Mean T2* error: %.2f%%\n\n', mean(T2_err));

    %% 10) Bar charts comparing GT vs LM
    quad_a_vec  = quad_a_lm(:)';
    quad_T2_vec = quad_T2_lm(:)';
    
    figure('Name','LM Fitting: Quadrant Comparison','Position',[100 100 1200 400]);
    % Amplitude comparison
    subplot(1,2,1);
    bar([s0; quad_a_vec]');
    xlabel('Quadrant');
    ylabel('a Value');
    title('True vs Estimated a');
    legend('True a','LM a','Location','best');
    grid on;

    % T2* comparison
    subplot(1,2,2);
    bar([T2; quad_T2_vec]');
    xlabel('Quadrant');
    ylabel('T_2^* (ms)');
    title('True vs Estimated T_2^*');
    legend('True T_2^*','LM T_2^*','Location','best');
    grid on;

    fprintf('Task 11 Synthetic LM Analysis Complete.\n');
end
