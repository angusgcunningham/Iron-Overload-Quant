function task_11_LM_model(patientPaths, consoleEstimates, roiRadius)
    % Task 11: Levenberg‑Marquardt vs ADMM vs Console T2* estimates
    % Uses four fixed ROIs (≈200 pixels) around (ref_row,ref_col) with offset.

    % ADMM values (from previous runs)
    admmEstimates = [10.29, 5.04, 15.15, 11.80];

    % --- ROI definition from Task 9 ---
    ref_row    = 120;                 % same as in Task 9
    ref_col    = 80;                  % same as in Task 9
    offset     = 8;                   % same as in Task 9
    roi_centers = [ ...
       ref_row-offset, ref_col-offset;   % top-left
       ref_row-offset, ref_col+offset;   % top-right
       ref_row+offset, ref_col-offset;   % bottom-left
       ref_row+offset, ref_col+offset ]; % bottom-right
    roi_radius = sqrt(200/pi);        % area ≈200 pixels

    % Preallocate results struct
    results = repmat(struct(...
        'PatientNumber',[],...
        'LM_T2star',[],...
        'ADMM_T2star',[],...
        'Console_T2star',[],...
        'Mean_a',zeros(1,4),...
        'Mean_T2',zeros(1,4)),...
        1, numel(patientPaths));

    % Create figure for amplitude/T2* for each patient
    figure('Name','LM Fitting Results: a & T2* Maps','Position',[50 50 1200 900]);

    for p = 1:numel(patientPaths)
        folder = patientPaths{p};
        fprintf('Processing patient %d: %s\n', p, folder);

        % --- Load DICOM series ---
        files = dir(fullfile(folder,'*.IMA'));
        if isempty(files)
            warning('No DICOM in %s, skipping...', folder);
            continue;
        end
        % read first echo to get dims
        firstFile = fullfile(folder, files(1).name);
        I0 = double(dicomread(firstFile));
        [rows, cols] = size(I0);
        TE = zeros(numel(files),1);
        data = zeros(rows,cols,numel(files));
        for i=1:numel(files)
            fn = fullfile(folder,files(i).name);
            info = dicominfo(fn);
            TE(i) = info.EchoTime;
            data(:,:,i) = double(dicomread(fn));
        end
        % sort by TE
        [TE,si] = sort(TE); data = data(:,:,si);

        % --- LM fit per pixel ---
        Npix = rows*cols;
        Y = reshape(data, rows*cols, [])';
        a_lm = nan(Npix,1);
        r_lm = nan(Npix,1);
        lb = [0;0]; ub = [Inf;0.5];
        opts = optimset('Display','off');
        model = @(p,t) p(1)*exp(-p(2)*t);
        for pix=1:Npix
            sig = Y(:,pix);
            if max(sig)<30, continue; end  % background
            init = [sig(1),0.05];
            pp = lsqcurvefit(model,init,TE,sig,lb,ub,opts);
            a_lm(pix)=pp(1);
            r_lm(pix)=pp(2);
        end

        % reshape & clamp
        a_img = reshape(a_lm, rows, cols);
        r_img = reshape(r_lm, rows, cols);
        r_img(r_img<eps)=eps;  % enforce positivity
        T2_img = 1./r_img;
        % median filter
        a_img = medfilt2(a_img,[3 3]);
        T2_img = medfilt2(T2_img,[3 3]);

        % --- ROI analysis: four fixed circles ---
        [X,Ygrid] = meshgrid(1:cols,1:rows);
        mean_a = nan(1,4);
        mean_T2 = nan(1,4);
        for k=1:4
            cr = roi_centers(k,1);
            cc = roi_centers(k,2);
            mask = ((Ygrid-cr).^2 + (X-cc).^2) <= roi_radius^2;
            mean_a(k)  = mean(a_img(mask));
            mean_T2(k) = mean(T2_img(mask & T2_img>0 & T2_img<100)); % ignore extreme
        end

        % store results
        results(p).PatientNumber  = p;
        results(p).LM_T2star      = mean_T2;
        results(p).ADMM_T2star    = admmEstimates(p);
        results(p).Console_T2star = consoleEstimates(p);
        results(p).Mean_a         = mean_a;
        results(p).Mean_T2        = mean_T2;

        % --- Plot amplitude + T2* side by side ---
        % Amplitude
        subplot(numel(patientPaths),2,(p-1)*2+1);
        imagesc(a_img,prctile(a_img(a_img>0),[1 99]));
        colormap(gca,hot); colorbar; axis image off;
        title(sprintf('Estimated a for Patient %d',p));
        hold on;
        for k=1:4
            theta=linspace(0,2*pi,100);
            plot(roi_centers(k,2)+roi_radius*cos(theta), roi_centers(k,1)+roi_radius*sin(theta),'w','LineWidth',1);
        end
        hold off;

        % T2*
        subplot(numel(patientPaths),2,(p-1)*2+2);
        T2c = T2_img;
        T2c(T2c>50)=50; T2c(T2c<0)=NaN;
        imagesc(T2c,[0 50]);
        colormap(gca,hot); colorbar; axis image off;
        title(sprintf('Patient %d - T2* Map (LM) mean: %.1f ms',p,mean_T2(1)));
        hold on;
        for k=1:4
            plot(roi_centers(k,2)+roi_radius*cos(theta), roi_centers(k,1)+roi_radius*sin(theta),'w','LineWidth',1);
        end
        hold off;
    end

    % --- Summary table ---
    fprintf('\nPatient |   ROI1   ROI2   ROI3   ROI4  |  ADMM  |  Console\n');
    for p=1:numel(results)
        mLM   = results(p).Mean_T2;
        mADMM = results(p).ADMM_T2star;
        mC    = results(p).Console_T2star;
        fprintf('  %2d    | %6.2f %6.2f %6.2f %6.2f | %6.2f | %6.2f\n',...
                p, mLM(1),mLM(2),mLM(3),mLM(4), mADMM, mC);
    end

    % --- Comparison bar chart ---
    figure('Name','T2* Comparison','Position',[200 200 600 300]);
    LMvals  = cell2mat({results(p).Mean_T2}');
    ADMMval = results(1).ADMM_T2star; % note ADMM is scalar per patient
    CONSval = results(1).Console_T2star;
    % build [patients x methods]
    comp = [mean(LMvals,2), ADMMval', CONSval'];
    bar(comp);
    set(gca,'XTick',1:numel(patientPaths),'XTickLabel',patientPaths);
    legend('LM','ADMM','Console','Location','best');
    ylabel('Mean T2* (ms)');
    title('Mean T2* in 4 ROIs');
end
